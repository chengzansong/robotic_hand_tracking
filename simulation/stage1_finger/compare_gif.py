"""Side-by-side GIFs of two designs over a pressure sweep, each frame annotated.

    python -m simulation.stage1_finger.compare_gif \\
        --free  simulation/outputs/eco30_sweep_10kPa/<run>  simulation/outputs/cloth_hi_sweep_meshB/<run> \\
        --blocked simulation/outputs/sw6_blocked_eco30_meshB/<run> simulation/outputs/cloth_hi_blocked_meshB/<run> \\
        --labels "Baseline: 4 mm Ecoflex 00-30 base" "Cloth layer in the base" \\
        --max-kPa 5 --out simulation/outputs/compare_cloth

For each pressure level solved in BOTH free-bending runs (0 included, up to
--max-kPa), renders each design in the same fixed side view, cut on the
centre plane so the chambers and any cloth layer show (the cloth part is drawn
in orange, everything else colored by max principal stretch on one shared
range), and writes:

    <out>/a.gif, <out>/b.gif     one GIF per design, identical camera and frames
    <out>/side_by_side.gif       the two frames next to each other
    <out>/frames/                the PNG frames
    <out>/stats.csv              the annotated numbers

Each frame shows: pressure; bend angle (tip-face rotation, metrics.py); turn
radius R and curvature 1000 / R in 1/m of the base underside (curvature.py:
R_fit when the fit is reliable, past about 90 deg, else the mean radius R_avg);
and the blocked tip force at that pressure from the matching blocked-force run,
linearly interpolated between its solved levels (and through 0 at 0 kPa), with
"interp." when the pressure is not a solved blocked level. The blocked run's
mesh name is shown when it differs from the free run's.

Frames are solved levels only; nothing is interpolated in shape.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .analysis.curvature import turn_radius_sweep
from .analysis.render import (BAR, _far_half, _gif, _mid_y, _pv, _side, _skin, _text_box,
                              _xz_bounds, volume_grid)
from .mesh.fe_mesh import FEMesh, mirror_half, principal_stretches

LAYER_COLOR = "#e6550d"
LAYER_PARTS = ("cloth",)
SIZE = (1200, 720)


def _mesh_name(run_dir):
    try:
        cfg = json.loads((Path(run_dir) / "metadata.json").read_text())["study_config"]["mesh"]
    except (OSError, KeyError, ValueError):
        return "?"
    g, w = cfg.get("size_mm"), cfg.get("wall_size_mm")
    return {(4.0, 1.5): "mesh A", (3.0, 1.0): "mesh B", (2.5, 0.8): "mesh C"}.get(
        (g, w), f"mesh {g:g}/{w:g} mm")


def load_free(run_dir):
    """Half mesh, pressures, displacements and per-level stats of a free run."""
    run_dir = Path(run_dir)
    mesh = FEMesh.load(run_dir / "mesh" / "mesh.npz")
    d = np.load(run_dir / "displacements.npz")
    p, u = d["pressure_kPa"], d["u"]
    res = json.loads((run_dir / "results.json").read_text())
    bend = {round(s["pressure_kPa"], 6): s["bend_angle_deg"] for s in res["steps"]}
    turn, _ = turn_radius_sweep(mesh, u, p)
    tr = {round(t["pressure_kPa"], 6): t for t in turn}
    stats = []
    for pk in p:
        key = round(float(pk), 6)
        t = tr.get(key)
        R = None
        if t:
            R = t["R_fit_mm"] if t["fit_reliable"] else t["R_avg_mm"]
        stats.append({"pressure_kPa": float(pk), "bend_deg": bend.get(key, 0.0),
                      "R_mm": R, "R_kind": ("fit" if t and t["fit_reliable"] else "avg") if t else "",
                      "kappa_per_m": 1000.0 / R if R else 0.0,
                      "underside_length": t["arc_stretch"] if t else 1.0})
    return {"dir": run_dir, "mesh": mesh, "p": np.asarray(p, float), "u": u, "stats": stats,
            "mesh_name": _mesh_name(run_dir)}


def load_blocked(run_dir):
    """(pressures, blocked forces in N, mesh name) of a blocked-force run, 0 included."""
    if not run_dir:
        return None
    res = json.loads((Path(run_dir) / "results.json").read_text())
    pts = sorted({(round(s["pressure_kPa"], 6), s.get("blocked_force_N") or 0.0)
                  for s in res["steps"]})
    p = np.array([a for a, _ in pts], float)
    f = np.array([b for _, b in pts], float)
    if p[0] > 0:
        p, f = np.r_[0.0, p], np.r_[0.0, f]
    return {"p": p, "F": f, "mesh_name": _mesh_name(run_dir),
            "reached": res.get("last_converged_pressure_kPa")}


def blocked_force_at(blocked, p_kPa):
    """(force N or None, exact solved level?) at a pressure; None past the last level."""
    if blocked is None or p_kPa > blocked["p"][-1] + 1e-9:
        return None, False
    exact = bool(np.any(np.isclose(blocked["p"], p_kPa, atol=1e-6)))
    return float(np.interp(p_kPa, blocked["p"], blocked["F"])), exact


def common_levels(a, b, max_kPa=None):
    pa = {round(float(x), 6) for x in a["p"]}
    pb = {round(float(x), 6) for x in b["p"]}
    lv = sorted(pa & pb)
    return [x for x in lv if max_kPa is None or x <= max_kPa + 1e-9]


def _index(run, p):
    return int(np.argmin(np.abs(run["p"] - p)))


def annotate(label, s, force, exact, blocked, free_mesh):
    lines = [label, f"p = {s['pressure_kPa']:.3g} kPa", f"bend {s['bend_deg']:.0f} deg"]
    if s["R_mm"] and s["bend_deg"] > 5:
        lines.append(f"turn radius {s['R_mm']:.1f} mm ({s['R_kind']}), "
                     f"curvature {s['kappa_per_m']:.1f} 1/m")
    else:
        lines.append("turn radius -, curvature 0 1/m")
    if force is None:
        lines.append("blocked tip force: not run" if blocked is None
                     else f"blocked tip force: not solved past {blocked['p'][-1]:g} kPa")
    else:
        tag = "" if exact else " (interp.)"
        mesh_tag = f", {blocked['mesh_name']}" if blocked["mesh_name"] != free_mesh else ""
        lines.append(f"blocked tip force {1000 * force:.1f} mN{tag}{mesh_tag}")
    return "\n".join(lines)


def _frame(pv, mesh, pts, lam, clim, bounds, text, path):
    """Side view of the far half (y > centre plane), cut face toward the camera."""
    grid = volume_grid(mesh, pts)
    grid.cell_data["stretch"] = lam
    y_mid = _mid_y(mesh)
    names = list(mesh.parts)
    layer_ids = [i for i, n in enumerate(names) if n in LAYER_PARTS]
    body = grid.extract_cells(np.flatnonzero(~np.isin(grid.cell_data["part"], layer_ids)))
    pl = pv.Plotter(off_screen=True, window_size=SIZE)
    pl.set_background("white")
    pl.add_mesh(_skin(_far_half(body.cell_data_to_point_data(), y_mid)), scalars="stretch",
                clim=clim, cmap="viridis",
                scalar_bar_args={"title": "max principal stretch", **BAR})
    if layer_ids:
        layer = grid.extract_cells(np.flatnonzero(np.isin(grid.cell_data["part"], layer_ids)))
        pl.add_mesh(_skin(_far_half(layer, y_mid)), color=LAYER_COLOR)
    _text_box(pl.add_text(text, font_size=11, position="upper_left"))
    _side(pl, bounds, aspect=SIZE[0] / SIZE[1])
    pl.screenshot(str(path))
    pl.close()


def make_compare(free_a, free_b, out, blocked_a=None, blocked_b=None,
                 labels=("A", "B"), max_kPa=None, clim=None, frame_s=0.8, hold_s=2.5):
    """Write the GIFs, frames and stats.csv. Returns {name: path}."""
    pv = _pv()
    out = Path(out)
    (out / "frames").mkdir(parents=True, exist_ok=True)
    runs = [load_free(free_a), load_free(free_b)]
    blocked = [load_blocked(blocked_a), load_blocked(blocked_b)]
    levels = common_levels(runs[0], runs[1], max_kPa)
    if len(levels) < 2:
        raise ValueError(f"the two free runs share too few pressure levels: {levels}")
    full = []
    for r in runs:
        m, u = (mirror_half(r["mesh"], r["u"]) if r["mesh"].info.get("symmetry")
                else (r["mesh"], r["u"]))
        full.append((m, u))
    frames_pts, lams = [], []
    for (m, u), r in zip(full, runs):
        pts = [m.nodes + u[_index(r, p)] for p in levels]
        frames_pts.append(pts)
        lams.append([principal_stretches(m.nodes, x, m.all_tets())[:, -1] for x in pts])
    if clim is None:
        top = max(float(l.max()) for ls in lams for l in ls)
        clim = (1.0, max(1.1, np.ceil(top * 10) / 10))
    xmin, xmax, zmin, zmax = _xz_bounds(np.vstack([x for pts in frames_pts for x in pts]))
    bounds = (xmin, xmax, zmin, zmax + 0.3 * (zmax - zmin))    # headroom for the text box

    rows, pngs = [], {"a": [], "b": []}
    for k, p in enumerate(levels):
        for j, key in enumerate(("a", "b")):
            r, (m, _) = runs[j], full[j]
            s = r["stats"][_index(r, p)]
            force, exact = blocked_force_at(blocked[j], p)
            text = annotate(labels[j], s, force, exact, blocked[j], r["mesh_name"])
            path = out / "frames" / f"{key}_{k:02d}_{p:g}kPa.png"
            _frame(pv, m, frames_pts[j][k], lams[j][k], clim, bounds, text, path)
            pngs[key].append(path)
            rows.append({"design": labels[j], "run": str(r["dir"]), **s,
                         "blocked_force_mN": None if force is None else 1000 * force,
                         "blocked_exact_level": exact})
    written = {key: _gif(pngs[key], out / f"{key}.gif", frame_s, hold_s) for key in pngs}
    from PIL import Image
    pairs = []
    for k, (fa, fb) in enumerate(zip(pngs["a"], pngs["b"])):
        ia, ib = Image.open(fa).convert("RGB"), Image.open(fb).convert("RGB")
        both = Image.new("RGB", (ia.width + ib.width, max(ia.height, ib.height)), "white")
        both.paste(ia, (0, 0))
        both.paste(ib, (ia.width, 0))
        path = out / "frames" / f"pair_{k:02d}.png"
        both.save(path)
        pairs.append(path)
    written["side_by_side"] = _gif(pairs, out / "side_by_side.gif", frame_s, hold_s)
    with open(out / "stats.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    written["stats"] = out / "stats.csv"
    return written


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--free", nargs=2, required=True, metavar=("RUN_A", "RUN_B"),
                    help="free-bending run folders")
    ap.add_argument("--blocked", nargs=2, metavar=("RUN_A", "RUN_B"),
                    help="blocked-force run folders for the same designs ('-' for none)")
    ap.add_argument("--labels", nargs=2, default=("Baseline", "Cloth base"))
    ap.add_argument("--max-kPa", type=float, default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    blk = [None if (b in (None, "-")) else b for b in (a.blocked or (None, None))]
    for name, path in make_compare(a.free[0], a.free[1], a.out, blk[0], blk[1],
                                   labels=a.labels, max_kPa=a.max_kPa).items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
