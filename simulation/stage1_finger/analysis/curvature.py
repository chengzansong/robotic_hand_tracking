"""Turn radius of the bent finger, from the deformed base underside.

At higher pressures the base curls into a rough circle. This module measures
that circle on the object-facing surface, the base underside (z = z_min) along
the finger centreline (the half model's symmetry plane), in the bending plane
spanned by the finger axis a and the bend direction b.

Definitions, per pressure level, over the chambered span (the undeformed
underside minus `margin_mm` at each end, so the end walls are excluded):

- R_fit_mm: least-squares circle through the deformed underside points
  (algebraic Kasa fit, then Gauss-Newton on the geometric distance).
- circularity_rms_pct: RMS distance of the points from that circle, as % of
  R_fit. Small means the underside really is a circular arc.
- R_avg_mm: deformed arc length / tangent turning over the span, the exact
  mean radius of curvature (R = s / theta). It matches the targets note,
  which uses R = L / theta.
- span_bend_deg: tangent turning over the span (the tip-normal bend angle in
  metrics.py also includes the end walls, so it reads a little higher).
- kappa_profile: mean curvature (1/mm) in sliding windows `window_mm` long
  (undeformed), centred every `step_mm`. It shows where bending concentrates.
- arc_stretch: deformed / undeformed underside length over the span.
- fit_reliable: False below `min_bend_deg` of span turning, where a circle
  fit through a nearly straight curve is ill-conditioned.
"""

import numpy as np

from .metrics import axis_vector

MIN_BEND_DEG = 90.0


def underside_centerline(mesh, axis="x", bend_toward="-z", tol_mm=1e-3):
    """Node ids on the base underside along the centreline, sorted along the axis.

    Uses the symmetry plane of a half model. For a full model, falls back to the
    underside nodes closest to the mid-width line (within 0.75 mm).
    """
    a, b = axis_vector(axis), axis_vector(bend_toward)
    w = np.cross(a, b)                       # width direction
    p = mesh.nodes
    hb = p @ b                               # height measured toward the bend side
    under = np.flatnonzero(hb > hb.max() - tol_mm)
    if "symmetry" in mesh.node_sets:
        sym = np.asarray(mesh.node_sets["symmetry"])
        y_mid = float(np.median(p[sym] @ w))
        band = tol_mm
    else:
        y_mid = float(0.5 * ((p @ w).min() + (p @ w).max()))
        band = 0.75
    ids = under[np.abs(p[under] @ w - y_mid) <= band]
    if len(ids) < 5:
        raise ValueError("too few underside centreline nodes for a turn-radius fit")
    ids = ids[np.argsort(p[ids] @ a)]
    if band > tol_mm:                         # full model: one node per axial position
        _, keep = np.unique(np.round(p[ids] @ a, 3), return_index=True)
        ids = ids[np.sort(keep)]
    return ids


def fit_circle(pts):
    """(center (2,), R, rms) of the best-fit circle through 2D points."""
    x, z = pts[:, 0], pts[:, 1]
    A = np.column_stack([x, z, np.ones_like(x)])
    D, E, F = np.linalg.lstsq(A, -(x ** 2 + z ** 2), rcond=None)[0]
    c = np.array([-D / 2, -E / 2])
    R = float(np.sqrt(max(c @ c - F, 1e-12)))
    for _ in range(20):                       # Gauss-Newton on geometric distance
        d = pts - c
        r = np.linalg.norm(d, axis=1)
        res = r - R
        J = np.column_stack([-d[:, 0] / r, -d[:, 1] / r, -np.ones_like(r)])
        step = np.linalg.lstsq(J, -res, rcond=None)[0]
        c, R = c + step[:2], R + step[2]
        if np.abs(step).max() < 1e-9 * max(R, 1.0):
            break
    rms = float(np.sqrt(np.mean((np.linalg.norm(pts - c, axis=1) - R) ** 2)))
    return c, float(abs(R)), rms


def _turning(xy):
    """Tangent turning (rad) of a polyline and the arc length it turns over.

    Segment angles are unwrapped along the curve. Each end tangent averages the
    two end segments, which puts it at the second and second-to-last points,
    so the arc length returned runs between those two points.
    """
    seg = np.diff(xy, axis=0)
    ds = np.linalg.norm(seg, axis=1)
    if len(seg) < 4:
        return 0.0, float(ds.sum())
    ang = np.unwrap(np.arctan2(seg[:, 1], seg[:, 0]))
    return float(ang[-2:].mean() - ang[:2].mean()), float(ds[1:-1].sum())


def turn_radius(mesh, u, ids=None, axis="x", bend_toward="-z", margin_mm=3.0,
                window_mm=30.0, step_mm=10.0, min_bend_deg=MIN_BEND_DEG):
    """Turn-radius metrics for one displacement field u (N, 3) mm."""
    a, b = axis_vector(axis), axis_vector(bend_toward)
    if ids is None:
        ids = underside_centerline(mesh, axis, bend_toward)
    s0 = mesh.nodes[ids] @ a
    lo, hi = s0.min() + margin_mm, s0.max() - margin_mm
    sel = (s0 >= lo - 1e-9) & (s0 <= hi + 1e-9)
    x = mesh.nodes[ids] + u[ids]
    xy = np.column_stack([x @ a, x @ b])     # bending plane; +b is toward the palm
    span = xy[sel]
    turn, s_turn = _turning(span)
    s_def = float(np.linalg.norm(np.diff(span, axis=0), axis=1).sum())
    s_ref = float(s0[sel].max() - s0[sel].min())
    c, R, rms = fit_circle(span)
    profile = []
    for mid in np.arange(lo + window_mm / 2, hi - window_mm / 2 + 1e-9, step_mm):
        w = (s0 >= mid - window_mm / 2 - 1e-9) & (s0 <= mid + window_mm / 2 + 1e-9)
        tw, sw = _turning(xy[w])
        profile.append({"x_mm": float(mid), "kappa_per_mm": tw / sw if sw else 0.0})
    bend = float(np.degrees(turn))
    return {
        "R_fit_mm": R,
        "center_mm": c.tolist(),
        "circularity_rms_pct": 100 * rms / R if R else None,
        "R_avg_mm": s_turn / turn if abs(turn) > 1e-9 else None,
        "span_bend_deg": bend,
        "arc_length_mm": s_def,
        "arc_stretch": s_def / s_ref,
        "span_mm": [float(lo), float(hi)],
        "kappa_profile": profile,
        "fit_reliable": bool(bend >= min_bend_deg),
    }


def turn_radius_sweep(mesh, u_all, pressures_kPa, **kw):
    """Per-level turn-radius rows; levels at zero pressure are skipped."""
    ids = underside_centerline(mesh, kw.get("axis", "x"), kw.get("bend_toward", "-z"))
    rows = []
    for p, u in zip(pressures_kPa, u_all):
        if p <= 0:
            continue
        rows.append({"pressure_kPa": float(p), **turn_radius(mesh, u, ids=ids, **kw)})
    return rows, ids
