"""STEP assembly -> tagged quadratic tetrahedral mesh, with gmsh.

Steps:
1. Import the design's assembly STEP in mm and match each solid to a design
   part by volume (each part STEP is imported in its own frame for that).
2. Fragment the solids so touching parts share a conformal interface. Any
   overlap becomes its own pieces, which go to `fem_tags.overlap_owner`.
3. Identify surfaces with geometric selectors from `[fem_tags]` in
   design.toml, never raw gmsh tags:
   - cavity: exterior surfaces grouped into connected shells by shared
     curves; every shell except the outer one bounds an internal cavity.
     With no inner shell the cavity is open, and surfaces strictly inside the
     overall bounding box are used instead.
   - clamp / tip: exterior planar faces lying on the bounding-box plane at the
     proximal / distal end of the finger axis.
   - slot walls: exterior planar faces normal to the finger axis, away from
     the ends, paired when one faces +axis and the next one along faces
     -axis within `slot_max_gap_mm` (the gap between neighbouring chambers).
   - symmetry: with `MeshSpec.symmetry = "half_y"` the solids are cut at the
     bounding-box mid-plane in y before fragmenting; the cut faces are the
     symmetry surface, and they are left out of the shell search.
4. Mesh with a size field that refines toward the cavity walls, raise to
   second order, and export an FEMesh plus a quality report.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .fe_mesh import (FEMesh, GMSH_TO_VTK_TET10, cavity_volume, orient_out_of_solid,
                      tet10_jacobians, tet10_volumes)

AXES = {"x": 0, "y": 1, "z": 2}

DEFAULT_TAGS = {
    "axis": "x",                 # finger length axis
    "clamp_end": "min",          # proximal end on that axis
    "bend_toward": "-z",         # side the finger curls toward (the base)
    "plane_tol_mm": 0.01,
    "cavity": "inner_shells",
    "interior_margin_mm": 0.5,   # only for the open-cavity fallback
    "overlap_owner": "top",
    "slot_walls": "facing_planes",
    "slot_max_gap_mm": 5.0,
}


@dataclass
class MeshSpec:
    size_mm: float = 2.0             # global target edge length
    wall_size_mm: float = 0.8        # target at the cavity walls
    wall_distance_mm: float = 2.0    # distance over which the size ramps up
    order: int = 2
    optimize_netgen: bool = True     # raises the worst-element quality noticeably
    symmetry: str = ""               # "" (full model) or "half_y" (cut at the y mid-plane)
    wall_elems: float = 0.0          # > 0: wall size = min(wall_size_mm, thinnest wall / wall_elems)
    # Target size inside named parts, e.g. {"cloth": 1.0} for a thin layer; the
    # smaller of this and the wall-refinement field applies there. Empty = off.
    part_size_mm: dict = None


def thinnest_wall_mm(design):
    """Thinnest wall of a design: fem_tags.min_wall_mm, else the smallest of the
    generated [parametric] wall, side, roof, and end thicknesses, else None."""
    if design.fem_tags.get("min_wall_mm"):
        return float(design.fem_tags["min_wall_mm"])
    import tomllib
    with open(design.directory / "design.toml", "rb") as f:
        par = tomllib.load(f).get("parametric", {})
    vals = [par[k] for k in ("wall_mm", "side_mm", "roof_mm", "end_mm") if k in par]
    return float(min(vals)) if vals else None


def resolve_wall_size(design, spec):
    """Effective cavity-wall mesh size for a design and MeshSpec."""
    if spec.wall_elems > 0:
        t = thinnest_wall_mm(design)
        if t is None:
            raise ValueError(f"design {design.id}: wall_elems needs fem_tags.min_wall_mm "
                             f"or a [parametric] section")
        return min(spec.wall_size_mm, t / spec.wall_elems)
    return spec.wall_size_mm


def fem_tags(design):
    tags = dict(DEFAULT_TAGS)
    tags.update(design.fem_tags or {})
    return tags


def _gmsh():
    import gmsh
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
    return gmsh


def _part_volumes(gmsh, design):
    """Volume (mm^3) of each design part, read from its own STEP file."""
    out = {}
    for part in design.parts:
        gmsh.model.add(f"part_{part.name}")
        ents = gmsh.model.occ.importShapes(str(part.step))
        gmsh.model.occ.synchronize()
        out[part.name] = sum(gmsh.model.occ.getMass(d, t) for d, t in ents if d == 3)
        gmsh.model.remove()
    return out


def _match_parts(gmsh, vols, part_volumes, rel_tol=1e-3):
    owner = {}
    for _, tag in vols:
        mass = gmsh.model.occ.getMass(3, tag)
        hits = [n for n, v in part_volumes.items() if abs(mass - v) <= rel_tol * v]
        if len(hits) != 1:
            raise ValueError(f"assembly solid {tag} ({mass:.1f} mm^3) matches parts {hits}; "
                             f"part volumes are {part_volumes}")
        owner[tag] = hits[0]
    if sorted(owner.values()) != sorted(part_volumes):
        raise ValueError(f"assembly solids {owner} do not map one-to-one onto parts")
    return owner


def _shells(gmsh, surfaces):
    """Group surfaces into connected components through shared curves."""
    parent = {s: s for s in surfaces}

    def find(s):
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s

    by_curve = {}
    for s in surfaces:
        for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
            by_curve.setdefault(abs(c), []).append(s)
    for group in by_curve.values():
        for s in group[1:]:
            parent[find(s)] = find(group[0])
    shells = {}
    for s in surfaces:
        shells.setdefault(find(s), []).append(s)
    return list(shells.values())


def _bbox(gmsh, dim, tags):
    boxes = np.array([gmsh.model.getBoundingBox(dim, t) for t in tags])
    return np.concatenate([boxes[:, :3].min(0), boxes[:, 3:].max(0)])


def _plane_faces(gmsh, faces, ax, value, tol):
    out = []
    for s in faces:
        b = gmsh.model.getBoundingBox(2, s)
        if b[ax + 3] - b[ax] < tol and abs(b[ax] - value) < tol:
            out.append(s)
    return out


def _slot_walls(gmsh, outer, vols, ax, box, tags):
    """Pairs of facing planar walls across the gaps between chambers."""
    tol = tags["plane_tol_mm"]
    lo_walls, hi_walls = {}, {}
    for s in outer:
        b = np.array(gmsh.model.getBoundingBox(2, s))
        a = b[ax]
        if b[ax + 3] - a > tol or abs(a - box[ax]) < tol or abs(a - box[ax + 3]) < tol:
            continue
        probe = (b[:3] + b[3:]) / 2
        probe[ax] -= 10 * tol
        solid_below = any(gmsh.model.isInside(3, v, probe.tolist()) for v in vols)
        (lo_walls if solid_below else hi_walls).setdefault(round(float(a), 6), []).append(s)
    pairs, used = [], set()
    for a in sorted(lo_walls):
        ahead = [h for h in sorted(hi_walls)
                 if a < h <= a + tags["slot_max_gap_mm"] and h not in used]
        if ahead:
            used.add(ahead[0])
            pairs.append((a, ahead[0], lo_walls[a], hi_walls[ahead[0]]))
    return pairs


def _select(gmsh, tags, volume_owner, sym_y=None):
    """Apply the geometric selectors. Returns surface tag lists and findings."""
    if tags["cavity"] != "inner_shells":
        raise ValueError(f"unknown fem_tags.cavity selector {tags['cavity']!r}; "
                         f"only 'inner_shells' is implemented")
    ax = AXES[tags["axis"]]
    tol = tags["plane_tol_mm"]
    vols = list(volume_owner)
    box = _bbox(gmsh, 3, vols)
    exterior = sorted({abs(t) for _, t in gmsh.model.getBoundary(
        [(3, v) for v in vols], combined=True, oriented=False)})
    symmetry = _plane_faces(gmsh, exterior, AXES["y"], sym_y, tol) if sym_y is not None else []
    exterior = [s for s in exterior if s not in symmetry]

    shells = _shells(gmsh, exterior)
    diag = [np.ptp(_bbox(gmsh, 2, sh).reshape(2, 3), axis=0).sum() for sh in shells]
    outer = shells[int(np.argmax(diag))]
    inner = [sh for sh in shells if sh is not outer]

    findings = {"n_exterior_surfaces": len(exterior), "n_shells": len(shells)}
    if inner:
        cavity = sorted(s for sh in inner for s in sh)
        findings["cavity_selector"] = f"{len(inner)} inner boundary shell(s)"
    else:
        m = tags["interior_margin_mm"]
        cavity = []
        for s in exterior:
            b = np.array(gmsh.model.getBoundingBox(2, s))
            if np.all(b[:3] > box[:3] + m) and np.all(b[3:] < box[3:] - m):
                cavity.append(s)
        findings["cavity_selector"] = (f"no inner shell; {len(cavity)} surfaces strictly "
                                       f"inside the bounding box by {m} mm")

    lo, hi = box[ax], box[ax + 3]
    clamp_at, tip_at = (lo, hi) if tags["clamp_end"] == "min" else (hi, lo)
    clamp = _plane_faces(gmsh, outer, ax, clamp_at, tol)
    tip = _plane_faces(gmsh, outer, ax, tip_at, tol)
    if not clamp or not tip:
        raise ValueError(f"no planar end faces found on {tags['axis']} = {clamp_at} or {tip_at}")

    interface = []
    for _, s in gmsh.model.getEntities(2):
        up, _ = gmsh.model.getAdjacencies(2, s)
        owners = {volume_owner[v] for v in up if v in volume_owner}
        if len(up) == 2 and len(owners) == 2:
            interface.append(s)

    groups = {"cavity": cavity, "clamp": clamp, "tip": tip, "interface": interface}
    slots = []
    if tags["slot_walls"] == "facing_planes":
        for k, (a, b, lo_s, hi_s) in enumerate(_slot_walls(gmsh, outer, vols, ax, box, tags), 1):
            name = f"slot{k:02d}"
            groups[f"{name}_lo"], groups[f"{name}_hi"] = lo_s, hi_s
            slots.append({"name": name, "lo_mm": a, "hi_mm": b, "gap_mm": round(b - a, 6)})
    elif tags["slot_walls"]:
        raise ValueError(f"unknown fem_tags.slot_walls selector {tags['slot_walls']!r}")
    if symmetry:
        groups["symmetry"] = symmetry

    findings.update(bbox_mm=[round(float(v), 4) for v in box],
                    clamp_plane=f"{tags['axis']} = {round(clamp_at, 4) + 0.0:g}",
                    tip_plane=f"{tags['axis']} = {round(tip_at, 4) + 0.0:g}",
                    slots=slots)
    return groups, findings


def _half_y(gmsh, vols, owner_of):
    """Keep y <= mid-plane of every solid. Returns new vols, owners, and the plane."""
    b = _bbox(gmsh, 3, [t for _, t in vols])
    y_mid = (b[1] + b[4]) / 2
    pad = 1.0
    new_vols, new_owner = [], {}
    for dim, tag in vols:
        box = gmsh.model.occ.addBox(b[0] - pad, b[1] - pad, b[2] - pad, b[3] - b[0] + 2 * pad,
                                    y_mid - b[1] + pad, b[5] - b[2] + 2 * pad)
        out, _ = gmsh.model.occ.intersect([(dim, tag)], [(3, box)])
        for d, t in out:
            new_vols.append((d, t))
            new_owner[t] = owner_of[tag]
    gmsh.model.occ.synchronize()
    return new_vols, new_owner, float(y_mid)


def _fragment(gmsh, vols, owner_of, overlap_owner):
    """Fragment all solids; return new volume -> part and an overlap report."""
    out, children = gmsh.model.occ.fragment(vols[:1], vols[1:])
    gmsh.model.occ.synchronize()
    parents = {}
    for (_, src), kids in zip(vols, children):
        for _, k in kids:
            parents.setdefault(k, []).append(owner_of[src])
    owner, overlap = {}, []
    for _, k in out:
        names = parents.get(k, [])
        if len(names) > 1:
            overlap.append(k)
            owner[k] = overlap_owner
        else:
            owner[k] = names[0]
    report = {"n_pieces": len(overlap),
              "volume_mm3": round(sum(gmsh.model.occ.getMass(3, k) for k in overlap), 4)}
    if overlap:
        b = _bbox(gmsh, 3, overlap)
        report["bbox_mm"] = [round(float(v), 4) for v in b]
        report["assigned_to"] = overlap_owner
    return owner, report


def _size_field(gmsh, cavity, spec, wall_size_mm, part_vols=None):
    """Wall-refinement threshold field, plus a constant size inside the parts
    of `part_vols` ({size_mm: [volume tags]}), combined by Min."""
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", cavity)
    gmsh.model.mesh.field.setNumber(f_dist, "Sampling", 20)
    f_th = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_th, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_th, "SizeMin", wall_size_mm)
    gmsh.model.mesh.field.setNumber(f_th, "SizeMax", spec.size_mm)
    gmsh.model.mesh.field.setNumber(f_th, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(f_th, "DistMax", spec.wall_distance_mm)
    fields = [f_th]
    for size, vols in (part_vols or {}).items():
        f_c = gmsh.model.mesh.field.add("Constant")
        gmsh.model.mesh.field.setNumber(f_c, "VIn", size)
        gmsh.model.mesh.field.setNumber(f_c, "VOut", spec.size_mm)
        gmsh.model.mesh.field.setNumbers(f_c, "VolumesList", vols)
        gmsh.model.mesh.field.setNumber(f_c, "IncludeBoundary", 1)
        fields.append(f_c)
    if len(fields) > 1:
        f_min = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(f_min, "FieldsList", fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(f_min)
    else:
        gmsh.model.mesh.field.setAsBackgroundMesh(f_th)
    for opt in ("MeshSizeFromPoints", "MeshSizeFromCurvature", "MeshSizeExtendFromBoundary"):
        gmsh.option.setNumber(f"Mesh.{opt}", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", spec.size_mm)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1 if spec.optimize_netgen else 0)


def _extract(gmsh, volume_owner, surf_groups):
    tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = coords.reshape(-1, 3)
    index = np.full(int(tags.max()) + 1, -1)
    index[tags.astype(int)] = np.arange(len(tags))

    parts, elem_tags = {}, []
    for vol, part in volume_owner.items():
        types, etags, conn = gmsh.model.mesh.getElements(3, vol)
        for ty, et, cn in zip(types, etags, conn):
            if ty != 11:
                raise ValueError(f"expected tet10 (gmsh type 11), got type {ty}")
            tets = index[cn.astype(int)].reshape(-1, 10)[:, GMSH_TO_VTK_TET10]
            parts.setdefault(part, []).append(tets)
            elem_tags.append(et)
    parts = {k: np.vstack(v) for k, v in parts.items()}

    surfaces = {}
    for name, ents in surf_groups.items():
        tris = []
        for s in ents:
            types, _, conn = gmsh.model.mesh.getElements(2, s)
            for ty, cn in zip(types, conn):
                if ty == 9:
                    tris.append(index[cn.astype(int)].reshape(-1, 6))
        surfaces[name] = np.vstack(tris) if tris else np.zeros((0, 6), int)

    used = np.unique(np.concatenate([t.ravel() for t in parts.values()]))
    if len(used) != len(coords):   # compact away geometry-only nodes
        remap = np.full(len(coords), -1)
        remap[used] = np.arange(len(used))
        coords = coords[used]
        parts = {k: remap[v] for k, v in parts.items()}
        surfaces = {k: remap[v] for k, v in surfaces.items()}
    return coords, parts, surfaces, np.concatenate(elem_tags)


def _quality(gmsh, elem_tags, nodes, parts):
    q = np.asarray(gmsh.model.mesh.getElementQualities(elem_tags, "minSICN"))
    gamma = np.asarray(gmsh.model.mesh.getElementQualities(elem_tags, "gamma"))
    tets = np.vstack(list(parts.values()))
    vol = tet10_volumes(nodes, tets)
    jac = tet10_jacobians(nodes, tets)
    return {
        "minSICN_min": round(float(q.min()), 4),
        "minSICN_p1": round(float(np.percentile(q, 1)), 4),
        "minSICN_median": round(float(np.median(q)), 4),
        "n_minSICN_below_0.2": int((q < 0.2).sum()),
        "gamma_min": round(float(gamma.min()), 4),
        "element_volume_min_mm3": float(vol.min()),
        "n_nonpositive_jacobian": int((jac <= 0).any(axis=1).sum()),
    }


def build_mesh(design, spec=MeshSpec(), out_dir=None):
    """Mesh a design's assembly. Returns (FEMesh, report dict)."""
    tags = fem_tags(design)
    if design.assembly is None:
        raise ValueError(f"design {design.id} has no assembly STEP; the mesher positions "
                         f"parts from the assembly")
    gmsh = _gmsh()
    try:
        part_volumes = _part_volumes(gmsh, design)
        gmsh.model.add(design.id)
        vols = [e for e in gmsh.model.occ.importShapes(str(design.assembly)) if e[0] == 3]
        gmsh.model.occ.synchronize()
        owner_of = _match_parts(gmsh, vols, part_volumes)
        raw_box = {owner_of[t]: [round(v, 4) for v in gmsh.model.getBoundingBox(3, t)]
                   for _, t in vols}

        sym_y = None
        if spec.symmetry == "half_y":
            vols, owner_of, sym_y = _half_y(gmsh, vols, owner_of)
        elif spec.symmetry:
            raise ValueError(f"unknown mesh symmetry {spec.symmetry!r}; use '' or 'half_y'")
        volume_owner, overlap = _fragment(gmsh, vols, owner_of, tags["overlap_owner"])
        groups, findings = _select(gmsh, tags, volume_owner, sym_y)

        for part in part_volumes:
            gmsh.model.addPhysicalGroup(3, [v for v, p in volume_owner.items() if p == part],
                                        name=part)
        for name, ents in groups.items():
            if ents:
                gmsh.model.addPhysicalGroup(2, ents, name=name)

        wall_size = resolve_wall_size(design, spec)
        part_vols = {}
        for part, size in (spec.part_size_mm or {}).items():
            if part not in part_volumes:
                raise ValueError(f"mesh part_size_mm names part {part!r}; design {design.id} "
                                 f"has {sorted(part_volumes)}")
            part_vols.setdefault(float(size), []).extend(
                v for v, p in volume_owner.items() if p == part)
        _size_field(gmsh, groups["cavity"], spec, wall_size, part_vols)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(spec.order)
        if out_dir is not None:
            Path(out_dir).mkdir(parents=True, exist_ok=True)
            gmsh.write(str(Path(out_dir) / "mesh.msh"))

        nodes, parts, surfaces, elem_tags = _extract(gmsh, volume_owner, groups)
        quality = _quality(gmsh, elem_tags, nodes, parts)
    finally:
        gmsh.finalize()

    all_tets = np.vstack(list(parts.values()))
    surfaces = {k: orient_out_of_solid(v, nodes, all_tets) if k != "interface" else v
                for k, v in surfaces.items()}
    node_sets = {k: np.unique(surfaces[k]) for k in ("clamp", "tip", "cavity", "symmetry")
                 if k in surfaces}

    v_cav, n_open = cavity_volume(nodes, surfaces["cavity"])
    if sym_y is not None:
        v_cav *= 2   # the cut cavity is capped on the symmetry plane; report the full model
    findings.update(
        part_volumes_mm3={k: round(v, 3) for k, v in part_volumes.items()},
        assembly_bbox_mm=raw_box,
        overlap=overlap,
        cavity_closed=n_open == 0,
        cavity_open_loops=n_open,
        cavity_volume_mm3=round(float(v_cav), 3),
    )
    symmetry = {"type": "half_y", "plane_y_mm": sym_y} if sym_y is not None else None
    mesh = FEMesh(nodes=nodes, parts=parts, surfaces=surfaces, node_sets=node_sets,
                  info={"design": design.id, "fem_tags": tags, "geometry": findings,
                        "symmetry": symmetry,
                        "slots": [sl["name"] for sl in findings["slots"]]})

    stats = {
        "n_nodes": int(len(nodes)),
        "n_elements": int(len(all_tets)),
        "elements_per_part": {k: int(len(v)) for k, v in parts.items()},
        "n_facets": {k: int(len(v)) for k, v in surfaces.items()},
        "n_node_sets": {k: int(len(v)) for k, v in node_sets.items()},
        "solid_volume_mm3": round(float(tet10_volumes(nodes, all_tets).sum()), 3),
        "extents_mm": [round(float(v), 4) for v in np.ptp(nodes, axis=0)],
        "spec": spec.__dict__,
        "wall_size_used_mm": wall_size,
        "thinnest_wall_mm": thinnest_wall_mm(design),
        "quality": quality,
    }
    mesh.info["stats"] = stats
    if out_dir is not None:
        mesh.save(Path(out_dir) / "mesh.npz")
        mesh.to_meshio().write(Path(out_dir) / "mesh.vtu")
    return mesh, {"stats": stats, "geometry": findings}
