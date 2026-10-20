"""Post-solve contact checks, independent of the solver.

- slot_gaps: minimum signed gap between the two facing walls of each slot
  between chambers. For every node of one wall, the distance to the nearest
  facet of the other wall is measured along that facet's outward normal.
  Positive means open, negative means the walls interpenetrate. Only nodes
  whose foot point lies within about one facet of the nearest facet centroid
  count, so wall regions that do not face each other are ignored.
- curl_flags: the finger curling back onto itself (tip onto the base or the
  clamp) is not modelled. This flags it: bend past `max_bend_deg`, or any
  tip-face node within `near_mm` of the proximal half of the finger.
"""

import numpy as np
from scipy.spatial import cKDTree

from ..mesh.fe_mesh import tri6_to_linear


def _facets(points, tris6):
    t = points[tri6_to_linear(tris6)]
    n = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    size = np.linalg.norm(t[:, 1] - t[:, 0], axis=1)
    return t.mean(1), n / np.linalg.norm(n, axis=1, keepdims=True), size


def _signed_gaps(points, a_tris6, b_tris6):
    """(node ids of wall a, signed gap of each to wall b), facing nodes only."""
    c, n, size = _facets(points, b_tris6)
    ids = np.unique(a_tris6)
    p = points[ids]
    _, k = cKDTree(c).query(p)
    d = np.einsum("ij,ij->i", p - c[k], n[k])
    lateral = np.linalg.norm(p - c[k] - d[:, None] * n[k], axis=1)
    keep = lateral <= 1.5 * size[k] + 1e-9
    return ids[keep], d[keep]


def _slot_pairs(mesh):
    for slot in mesh.info.get("slots", []):
        lo, hi = mesh.surfaces.get(f"{slot}_lo"), mesh.surfaces.get(f"{slot}_hi")
        if lo is not None and hi is not None and len(lo) and len(hi):
            yield slot, lo, hi


def slot_gaps(mesh, points):
    """{slot: min signed gap in mm} for the deformed node positions `points`."""
    out = {}
    for slot, lo, hi in _slot_pairs(mesh):
        g = np.concatenate([_signed_gaps(points, lo, hi)[1], _signed_gaps(points, hi, lo)[1]])
        out[slot] = float(g.min()) if len(g) else None
    return out


def wall_gap_field(mesh, points):
    """Per-node signed gap on slot-wall nodes (NaN elsewhere), for rendering."""
    field = np.full(len(points), np.nan)
    for _, lo, hi in _slot_pairs(mesh):
        for a, b in ((lo, hi), (hi, lo)):
            ids, d = _signed_gaps(points, a, b)
            field[ids] = np.fmin(field[ids], d)
    return field


def curl_flags(mesh, points, bend_angle_deg, axis=0, max_bend_deg=300.0, near_mm=1.0):
    ref = mesh.nodes[:, axis]
    proximal = np.flatnonzero(ref < ref.min() + 0.5 * np.ptp(ref))
    tip = mesh.node_sets["tip"]
    d, _ = cKDTree(points[proximal]).query(points[tip])
    return {
        "bend_over_limit": bool(abs(bend_angle_deg) > max_bend_deg),
        "tip_to_proximal_min_mm": float(d.min()),
        "tip_near_proximal": bool(d.min() < near_mm),
    }
