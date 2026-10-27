"""Turn solver output (times, displacement fields, reactions) into results.json rows.

Pressure at each output time is p = max_kPa * t, because the writer ramps
the cavity pressure linearly over t in [0, 1]. A zero-pressure row is
prepended.

Half models (mesh symmetry "half_y") are mirrored to the full model before
any metric, so volumes, tip position, and stretches are full-model values.
Logged reactions come from the half and are doubled in x and z; y cancels.

Blocked force: FEBio 4.13 logs the nodal reaction R as the force the body
exerts on the support (checked 2026-10-20 on the single-tet10 tension test,
where the held x = 0 face logs +x). So blocked_force_N = R_tip . b, positive
when the finger pushes toward its bend side b. The same holds for the
blocked_line case: only the support-line nodes carry reactions there, and the
tip log sums the whole face.
"""

import numpy as np

from ..mesh.fe_mesh import mirror_half
from .contact import curl_flags, slot_gaps
from .metrics import axis_vector, step_metrics, unwrap_angles


def _full_reaction(r, half):
    if r is None:
        return None
    r = np.asarray(r, float)
    return (r * [2.0, 0.0, 2.0]).tolist() if half else r.tolist()


def sweep_rows(mesh, times, u, raw_steps, load, axis="x", bend_toward="-z",
               penetration_tol_mm=0.05):
    b = axis_vector(bend_toward)
    half = bool((mesh.info or {}).get("symmetry"))
    times = np.concatenate([[0.0], np.asarray(times, float)])
    u = (np.concatenate([np.zeros((1,) + u.shape[1:]), u]) if len(u)
         else np.zeros((1,) + mesh.nodes.shape))
    raw_steps = [{"time": 0.0, "tip_reaction_N": None, "clamp_reaction_N": None}] + list(raw_steps)
    full, u_full = mirror_half(mesh, u) if half else (mesh, u)
    levels = {round(f, 6) for f in [0.0] + load.fractions()}
    ax = "xyz".index(axis.lstrip("+-"))

    rows = []
    for t, ui, raw in zip(times, u_full, raw_steps):
        row = {"time": float(t), "pressure_kPa": float(load.max_kPa * t),
               "requested_level": round(float(t), 6) in levels}
        row.update(step_metrics(full, ui, axis=axis, bend_toward=bend_toward))
        x = full.nodes + ui
        gaps = slot_gaps(full, x)
        row["slot_min_gap_mm"] = gaps
        row["min_slot_gap_mm"] = min((g for g in gaps.values() if g is not None), default=None)
        row["penetration"] = (row["min_slot_gap_mm"] is not None
                              and row["min_slot_gap_mm"] < -penetration_tol_mm)
        row["curl"] = curl_flags(full, x, row["bend_angle_deg"], axis=ax)
        row["tip_reaction_N"] = _full_reaction(raw.get("tip_reaction_N"), half)
        row["clamp_reaction_N"] = _full_reaction(raw.get("clamp_reaction_N"), half)
        if load.blocked:
            r = row["tip_reaction_N"]
            row["blocked_force_N"] = float(np.asarray(r) @ b) if r is not None else 0.0
        rows.append(row)
    unwrap_angles(rows)
    for r in rows:      # re-check the bend limit on the unwrapped angle
        r["curl"]["bend_over_limit"] = bool(abs(r["bend_angle_deg"]) > 300.0)
    v0 = rows[0]["cavity_volume_mm3"]
    for r in rows:
        r["cavity_volume_change_mm3"] = r["cavity_volume_mm3"] - v0
    return rows, u, np.array([r["pressure_kPa"] for r in rows])


def summarize(rows, load, solver_status):
    """Run-level flags for results.json."""
    last = rows[-1]
    pen = [r for r in rows if r["penetration"]]
    curl = [r for r in rows if r["curl"]["bend_over_limit"] or r["curl"]["tip_near_proximal"]]
    gaps = [r["min_slot_gap_mm"] for r in rows if r["min_slot_gap_mm"] is not None]
    return {
        "solver_status": solver_status,
        "reached_full_load": bool(last["time"] > 1 - 1e-6),
        "last_converged_pressure_kPa": last["pressure_kPa"],
        "requested_max_kPa": load.max_kPa,
        "min_slot_gap_mm": min(gaps) if gaps else None,
        "penetration_warning": bool(pen),
        "first_penetration_kPa": pen[0]["pressure_kPa"] if pen else None,
        "curl_warning": bool(curl),
        "first_curl_kPa": curl[0]["pressure_kPa"] if curl else None,
    }
