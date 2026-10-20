"""Solver-independent tetrahedral mesh container and element math.

Connectivity is 0-based. Quadratic tetrahedra use the VTK / FEBio node
order: corners 0-3, then edge mid-nodes (0,1) (1,2) (2,0) (0,3) (1,3) (2,3).
gmsh's Tetrahedron10 swaps the last two, so the mesher permutes them on import.
Quadratic triangles (tri6) use corners 0-2, then edges (0,1) (1,2) (2,0),
which gmsh, VTK, and FEBio share.

Surface facets in `surfaces` are wound so the right-hand normal points out of
the solid. On the cavity that means the normal points into the air, which is
what FEBio's pressure load expects for a positive inflating pressure.
"""

from dataclasses import dataclass, field
import json
from pathlib import Path

import numpy as np

GMSH_TO_VTK_TET10 = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]

# Local tet10 faces, wound outward, as (corners | mid-nodes). Same table as
# FEBio's FEElement::GetFace for tet10 (FECore/FEElement.cpp, v4.13).
TET10_FACES = np.array([
    [0, 1, 3, 4, 8, 7],
    [1, 2, 3, 5, 9, 8],
    [2, 0, 3, 6, 7, 9],
    [2, 1, 0, 5, 4, 6],
])


@dataclass
class FEMesh:
    nodes: np.ndarray                       # (N, 3) mm
    parts: dict                             # part name -> (E, 10) int tet10
    surfaces: dict = field(default_factory=dict)   # name -> (F, 6) int tri6
    node_sets: dict = field(default_factory=dict)  # name -> (n,) int
    info: dict = field(default_factory=dict)       # geometry findings, JSON-safe

    @property
    def n_elements(self):
        return sum(len(e) for e in self.parts.values())

    def all_tets(self):
        return np.vstack(list(self.parts.values()))

    def save(self, path):
        path = Path(path)
        arrays = {"nodes": self.nodes}
        for k, v in self.parts.items():
            arrays[f"part__{k}"] = v
        for k, v in self.surfaces.items():
            arrays[f"surface__{k}"] = v
        for k, v in self.node_sets.items():
            arrays[f"nodeset__{k}"] = v
        arrays["info_json"] = np.array(json.dumps(self.info))
        np.savez_compressed(path, **arrays)

    @classmethod
    def load(cls, path):
        data = np.load(path)
        mesh = cls(nodes=data["nodes"], parts={})
        for key in data.files:
            if key.startswith("part__"):
                mesh.parts[key[6:]] = data[key]
            elif key.startswith("surface__"):
                mesh.surfaces[key[9:]] = data[key]
            elif key.startswith("nodeset__"):
                mesh.node_sets[key[9:]] = data[key]
        mesh.info = json.loads(str(data["info_json"]))
        return mesh

    def to_meshio(self, points=None):
        """meshio Mesh with a 'part' cell label (0, 1, ... in parts order)."""
        import meshio
        cells, labels = [], []
        for i, tets in enumerate(self.parts.values()):
            cells.append(("tetra10", tets))
            labels.append(np.full(len(tets), i))
        return meshio.Mesh(self.nodes if points is None else points, cells,
                           cell_data={"part": labels})


# --- element math -----------------------------------------------------------

def tet10_shape_grad(r, s, t):
    """dN/d(r,s,t) for tet10 in VTK order, shape (10, 3)."""
    L = np.array([1 - r - s - t, r, s, t])
    dL = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    g = np.zeros((10, 3))
    for a in range(4):
        g[a] = (4 * L[a] - 1) * dL[a]
    for a, (i, j) in enumerate([(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)], start=4):
        g[a] = 4 * (L[i] * dL[j] + L[j] * dL[i])
    return g


CORNER_POINTS = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0.25, 0.25, 0.25)]


def tet10_jacobians(nodes, tets, points=CORNER_POINTS):
    """det(dX/dxi) at the given natural points, shape (E, len(points))."""
    X = nodes[tets]                                   # (E, 10, 3)
    dets = []
    for p in points:
        g = tet10_shape_grad(*p)                      # (10, 3)
        J = np.einsum("eai,aj->eij", X, g)
        dets.append(np.linalg.det(J))
    return np.stack(dets, axis=1)


def tet10_volumes(nodes, tets):
    """Exact volumes of straight-sided tets from the corners (mm^3)."""
    X = nodes[tets[:, :4]]
    return np.einsum("ei,ei->e", X[:, 1] - X[:, 0],
                     np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 0])) / 6.0


def principal_stretches(ref_nodes, def_nodes, tets):
    """Principal stretches at each tet10 centroid, shape (E, 3), ascending."""
    g = tet10_shape_grad(0.25, 0.25, 0.25)
    X = ref_nodes[tets]
    x = def_nodes[tets]
    J0 = np.einsum("eai,aj->eij", X, g)               # dX/dxi
    dNdX = np.einsum("aj,eji->eai", g, np.linalg.inv(J0))
    F = np.einsum("eai,eaj->eij", x, dNdX)
    C = np.einsum("eki,ekj->eij", F, F)
    return np.sqrt(np.clip(np.linalg.eigvalsh(C), 0, None))


# --- surfaces ---------------------------------------------------------------

def tri6_flip(tris):
    """Reverse winding of tri6 facets."""
    return tris[:, [0, 2, 1, 5, 4, 3]]


def orient_out_of_solid(tris, nodes, tets):
    """Wind each tri6 so its normal points away from the tet it bounds."""
    lookup = {}
    for e, tet in enumerate(tets[:, :4]):
        for a, b, c in ((0, 1, 3), (1, 2, 3), (2, 0, 3), (0, 2, 1)):
            lookup.setdefault(frozenset((tet[a], tet[b], tet[c])), e)
    out = tris.copy()
    missing = 0
    for i, tri in enumerate(tris):
        e = lookup.get(frozenset(tri[:3]))
        if e is None:
            missing += 1
            continue
        p = nodes[tri[:3]]
        n = np.cross(p[1] - p[0], p[2] - p[0])
        if np.dot(n, p.mean(0) - nodes[tets[e, :4]].mean(0)) < 0:
            out[i] = tri[[0, 2, 1, 5, 4, 3]]
    if missing:
        raise ValueError(f"{missing} surface facets do not match any tet face")
    return out


def tri6_to_linear(tris):
    """Split each tri6 into its 4 corner/mid-node sub-triangles, same winding."""
    return np.vstack([tris[:, [0, 3, 5]], tris[:, [3, 1, 4]],
                      tris[:, [5, 4, 2]], tris[:, [3, 4, 5]]])


def boundary_loops(tris3):
    """Edge loops where a triangle surface is open (edges used once)."""
    count = {}
    for t in tris3:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            key = (min(a, b), max(a, b))
            count[key] = count.get(key, 0) + 1
    # keep the directed edges whose undirected key is used once
    nxt = {}
    for t in tris3:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            if count[(min(a, b), max(a, b))] == 1:
                nxt[a] = b
    loops = []
    while nxt:
        start, cur = next(iter(nxt.items()))
        loop = [start]
        del nxt[start]
        while cur != start:
            loop.append(cur)
            cur = nxt.pop(cur)
        loops.append(loop)
    return loops


def enclosed_volume(points, tris3, caps=None):
    """Volume enclosed by triangles whose normals point OUT of the region.

    `caps` are boundary loops (node index lists) closed with a fan from the
    loop centroid, for a cavity that opens through an inlet port.
    """
    p = points[tris3]
    v = np.einsum("ei,ei->", p[:, 0], np.cross(p[:, 1], p[:, 2])) / 6.0
    for loop in caps or []:
        # a cap facet (c, b, a) runs the shared edge opposite to the surface
        c = points[loop].mean(0)
        a = points[loop]
        b = points[np.roll(loop, -1)]
        v += np.einsum("ei,ei->", np.broadcast_to(c, a.shape), np.cross(b, a)) / 6.0
    return v


def cavity_volume(points, cavity_tris6):
    """Air volume of the cavity, mm^3, by the divergence theorem.

    Cavity facets point out of the solid, which is into the air, so the
    region-outward orientation is the reverse. An open cavity (inlet port)
    is capped with fans before integration.
    """
    tris3 = tri6_to_linear(cavity_tris6)[:, ::-1]
    loops = boundary_loops(tris3)
    return enclosed_volume(points, tris3, caps=loops), len(loops)


# --- symmetry ----------------------------------------------------------------

# Reflection flips orientation; swapping corners 1 and 2 (and the matching
# mid-nodes) makes the mirrored tet10 positive again.
_TET10_REFLECT = [0, 2, 1, 3, 6, 5, 4, 7, 9, 8]


def mirror_half(mesh, u=None):
    """Full-model FEMesh (and displacements) from a half_y symmetry mesh.

    Nodes on the plane are duplicated, not merged. Surfaces keep their
    out-of-solid winding. The cut cavity of each half is open on the plane,
    and the fan caps of the two halves cancel, so cavity_volume() of the
    mirrored surface is the full closed volume.
    """
    sym = (mesh.info or {}).get("symmetry")
    if not sym:
        return (mesh, u) if u is not None else mesh
    y0 = sym["plane_y_mm"]
    n = len(mesh.nodes)
    refl = mesh.nodes.copy()
    refl[:, 1] = 2 * y0 - refl[:, 1]
    full = FEMesh(
        nodes=np.vstack([mesh.nodes, refl]),
        parts={k: np.vstack([v, v[:, _TET10_REFLECT] + n]) for k, v in mesh.parts.items()},
        surfaces={k: np.vstack([v, tri6_flip(v) + n]) for k, v in mesh.surfaces.items()},
        node_sets={k: np.concatenate([v, v + n]) for k, v in mesh.node_sets.items()},
        info={**mesh.info, "symmetry": None, "mirrored_from_half": sym},
    )
    if u is None:
        return full
    u = np.asarray(u)
    ur = u.copy()
    ur[..., 1] *= -1
    return full, np.concatenate([u, ur], axis=-2)
