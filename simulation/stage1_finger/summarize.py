"""Collect every finished run into one results table (no re-solving).

    python -m simulation.stage1_finger.summarize                 # all runs under simulation/outputs
    python -m simulation.stage1_finger.summarize --latest        # newest run per study only
    python -m simulation.stage1_finger.summarize --compare cmp_eco30_meshA cmp_ds20base_meshA         --out simulation/outputs/compare.png                     # overlay studies (newest run each)

Writes simulation/outputs/summary.csv (one row per run and pressure level) and
simulation/outputs/summary.md (a run index plus a per-level table). Runs saved
before the turn-radius metric existed get it computed from displacements.npz.
The curated tables in simulation/docs/experiments.md are copied from summary.md.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .analysis.curvature import turn_radius_sweep
from .mesh.fe_mesh import FEMesh
from .run import OUTPUTS_DIR

MESH_NAMES = {(4.0, 1.5): "A", (3.0, 1.0): "B", (2.5, 0.8): "C"}

COLUMNS = ["study", "run", "design", "top", "base", "mesh", "case", "pressure_kPa",
           "bend_deg", "R_fit_mm", "R_avg_mm", "circularity_pct", "arc_stretch",
           "blocked_force_N", "dV_pct", "stretch_top", "stretch_base", "min_gap_mm"]


def _mesh_name(cfg):
    m = cfg.get("mesh", {})
    key = (float(m.get("size_mm", 0)), float(m.get("wall_size_mm", 0)))
    name = MESH_NAMES.get(key, f"{key[0]:g}/{key[1]:g}")
    return name + (" half" if m.get("symmetry") else "")


def _turn_rows(run_dir, results):
    if results.get("turn_radius"):
        return results["turn_radius"]
    if results.get("load_case") != "free" or not (run_dir / "displacements.npz").exists():
        return []
    try:
        mesh = FEMesh.load(run_dir / "mesh" / "mesh.npz")
        d = np.load(run_dir / "displacements.npz")
        return turn_radius_sweep(mesh, d["u"], d["pressure_kPa"])[0]
    except (OSError, ValueError, KeyError):
        return []


def run_rows(run_dir):
    """(index row, per-level rows) for one run folder, or None if it has no results."""
    run_dir = Path(run_dir)
    try:
        res = json.loads((run_dir / "results.json").read_text())
        meta = json.loads((run_dir / "metadata.json").read_text())
    except (OSError, ValueError):
        return None
    cfg = meta.get("study_config", {})
    mats = cfg.get("materials", {})
    turn = {round(t["pressure_kPa"], 6): t for t in _turn_rows(run_dir, res)}
    v0 = res["steps"][0]["cavity_volume_mm3"]
    base = {"study": res["study"], "run": run_dir.name, "design": res["design"],
            "top": mats.get("top", ""), "base": mats.get("base", ""), "mesh": _mesh_name(cfg),
            "case": res["load_case"]}
    levels = []
    for s in res["steps"][1:]:
        t = turn.get(round(s["pressure_kPa"], 6), {})
        st = s["principal_stretch"]
        levels.append({**base,
                       "pressure_kPa": s["pressure_kPa"],
                       "bend_deg": s["bend_angle_deg"],
                       "R_fit_mm": t.get("R_fit_mm") if t.get("fit_reliable") else None,
                       "R_avg_mm": t.get("R_avg_mm"),
                       "circularity_pct": t.get("circularity_rms_pct") if t.get("fit_reliable") else None,
                       "arc_stretch": t.get("arc_stretch"),
                       "blocked_force_N": s.get("blocked_force_N"),
                       "dV_pct": 100 * (s["cavity_volume_mm3"] / v0 - 1),
                       "stretch_top": st.get("top", {}).get("max"),
                       "stretch_base": st.get("base", {}).get("max"),
                       "min_gap_mm": s.get("min_slot_gap_mm")})
    solve = meta.get("solve", {})
    index = {**base, "status": res.get("solver_status", ""),
             "levels": (f"{res['last_converged_pressure_kPa']:g} of {res['requested_max_kPa']:g} kPa"
                        if "last_converged_pressure_kPa" in res
                        else f"{res['steps'][-1]['pressure_kPa']:g} kPa"),
             "wall_min": round(solve.get("wall_time_s", 0) / 60, 1) if solve.get("wall_time_s") else "",
             "path": str(run_dir.relative_to(OUTPUTS_DIR.parent.parent)).replace("\\", "/")}
    return index, levels


def _fmt(v, nd=3):
    if v is None or v == "":
        return "—"
    if isinstance(v, float):
        return f"{v:.{nd}g}" if abs(v) < 1000 else f"{v:.0f}"
    return str(v)


def collect(root=OUTPUTS_DIR, latest=False):
    runs = sorted(p.parent for p in Path(root).glob("*/*/results.json"))
    if latest:
        newest = {}
        for r in runs:
            newest[r.parent.name] = r
        runs = sorted(newest.values())
    out = [x for x in (run_rows(r) for r in runs) if x]
    return [i for i, _ in out], [lv for _, lvs in out for lv in lvs]


def write(index, levels, root=OUTPUTS_DIR):
    root = Path(root)
    with open(root / "summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(levels)
    icols = ["study", "design", "top", "base", "mesh", "case", "status", "levels", "wall_min", "path"]
    lines = ["# Run summary (generated by simulation/stage1_finger/summarize.py)", "",
             "## Runs", "", "| " + " | ".join(icols) + " |", "|" + " --- |" * len(icols)]
    lines += ["| " + " | ".join(_fmt(r[c]) for c in icols) + " |" for r in index]
    lcols = COLUMNS[:3] + COLUMNS[5:]
    lines += ["", "## Per level", "", "| " + " | ".join(lcols) + " |", "|" + " --- |" * len(lcols)]
    lines += ["| " + " | ".join(_fmt(r[c]) for c in lcols) + " |" for r in levels]
    (root / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return root / "summary.csv", root / "summary.md"


def plot_compare(studies, out_path, root=OUTPUTS_DIR):
    """Bend, turn radius and blocked force vs pressure for the newest run of each study."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .analysis.plots import REFERENCE_OBJECTS

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    for name in studies:
        runs = sorted(Path(root).glob(f"{name}/*/results.json"))
        if not runs:
            print(f"no results for {name}")
            continue
        _, levels = run_rows(runs[-1].parent)
        p = [r["pressure_kPa"] for r in levels]
        if levels and levels[0]["case"] == "blocked_force":
            axs[2].plot(p, [1000 * r["blocked_force_N"] for r in levels], "o-", label=name)
            continue
        axs[0].plot(p, [r["bend_deg"] for r in levels], "o-", label=name)
        pr = [(r["pressure_kPa"], r["R_fit_mm"]) for r in levels if r["R_fit_mm"]]
        if pr:
            axs[1].plot(*zip(*pr), "o-", label=name)
    for label, _, size in REFERENCE_OBJECTS:
        axs[1].axhline(size / 2, ls="--", lw=0.8, color="0.5")
        axs[1].annotate(label, (axs[1].get_xlim()[0], size / 2), fontsize=7, color="0.4", va="bottom")
    axs[0].axhline(200, ls="--", lw=0.8, color="0.5")
    axs[0].set(xlabel="chamber pressure (kPa)", ylabel="bend angle (deg)", title="free bend (200° target dashed)")
    axs[1].set(xlabel="chamber pressure (kPa)", ylabel="turn radius R fit (mm)",
               title="turn radius (bends ≥ 90° only)")
    axs[2].set(xlabel="chamber pressure (kPa)", ylabel="blocked tip force (mN)", title="blocked tip force")
    for ax in axs:
        ax.grid(alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7)
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", default=str(OUTPUTS_DIR))
    ap.add_argument("--latest", action="store_true", help="newest run per study only")
    ap.add_argument("--compare", nargs="+", metavar="STUDY", help="overlay plot of these studies")
    ap.add_argument("--out", default=str(OUTPUTS_DIR / "compare.png"), help="path for --compare")
    a = ap.parse_args(argv)
    if a.compare:
        print(plot_compare(a.compare, a.out, a.root))
        return
    index, levels = collect(a.root, a.latest)
    for p in write(index, levels, a.root):
        print(p)


if __name__ == "__main__":
    main()
