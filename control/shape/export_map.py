"""Export a Stage 1 run's results.json as a committed shape-map table.

    python -m control.shape.export_map <results.json> -o control/shape/maps/<name>.json

Reads the run with json and numpy only; nothing from simulation/ is imported,
so control can use the table without the simulation code or a solver.

Per pressure level the table holds bend angle, curvature, cavity volume and
blocked force (null where the run does not have them). Curvature comes from the
run's `turn_radius` rows when it has them: 1/R_fit once the span turns past
about 90 degrees, 1/R_avg below that (R_fit is ill-conditioned for a nearly
straight finger). A level without a turn-radius row falls back to
bend angle / arc length. `curvature_source` records which rule gave each point.

Some older runs predate `turn_radius`. For those, `--radius-mm 2:58.4,3:37.1`
supplies R_fit values read from that run's turn-radius plot, and the arc
length for the remaining levels is taken as the median of bend * R over them.
"""

import argparse
import json
from pathlib import Path

import numpy as np

DEFAULT_ARC_MM = 106.0     # chambered span of finger_poc (x = 3..109 mm), used if nothing better


def parse_radii(text):
    out = {}
    for item in text.split(","):
        p, r = item.split(":")
        out[float(p)] = float(r)
    return out


def export_map(results_path, design=None, material=None, radius_mm=None, arc_mm=None):
    """Build the table dict from a results.json path."""
    results_path = Path(results_path)
    res = json.loads(results_path.read_text())
    steps = [s for s in res["steps"] if s.get("requested_level", True)]
    steps.sort(key=lambda s: s["pressure_kPa"])
    p = np.array([s["pressure_kPa"] for s in steps], float)
    bend = np.array([s["bend_angle_deg"] for s in steps], float)

    rows = {round(r["pressure_kPa"], 6): r for r in res.get("turn_radius") or []}
    radii = {round(k, 6): v for k, v in (radius_mm or {}).items()}
    if arc_mm is None:
        implied = [np.radians(b) * radii[round(pp, 6)] for pp, b in zip(p, bend)
                   if round(pp, 6) in radii]
        arc_mm = float(np.median(implied)) if implied else DEFAULT_ARC_MM

    kappa, source = [], []
    for pp, b in zip(p, bend):
        key = round(float(pp), 6)
        row = rows.get(key)
        if pp <= 0:
            kappa.append(0.0)
            source.append("zero")
        elif row and row["R_fit_mm"] and row.get("fit_reliable", b >= 90):
            kappa.append(1000.0 / row["R_fit_mm"])
            source.append("1/R_fit")
        elif row and row["R_avg_mm"]:
            kappa.append(1000.0 / row["R_avg_mm"])
            source.append("1/R_avg")
        elif key in radii:
            kappa.append(1000.0 / radii[key])
            source.append("1/R_fit (from --radius-mm)")
        else:
            kappa.append(float(np.radians(b)) / (arc_mm / 1000.0))
            source.append(f"bend/arc_length ({arc_mm:.1f} mm)")

    vol = [s.get("cavity_volume_mm3") for s in steps]
    force = [s.get("blocked_force_N") for s in steps]
    return {
        "design": design or res.get("design"),
        "material": material,
        "source_run": results_path.as_posix(),
        "study": res.get("study"),
        "load_case": res.get("load_case"),
        "solver_status": res.get("solver_status"),
        "curvature_rule": "1/R_fit past ~90 deg of span turning, 1/R_avg below, else bend/arc_length; see curvature_source",
        "curvature_source": source,
        "arc_length_mm": arc_mm,
        "pressure_kPa": p.tolist(),
        "bend_deg": bend.tolist(),
        "curvature_per_m": [round(k, 4) for k in kappa],
        "cavity_volume_mm3": None if any(v is None for v in vol) else vol,
        "blocked_force_N": None if all(f is None for f in force) else force,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("results", help="path to a Stage 1 results.json")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--design")
    ap.add_argument("--material")
    ap.add_argument("--radius-mm", type=parse_radii,
                    help="kPa:R_fit_mm pairs for runs without turn_radius, e.g. 2:58.4,3:37.1")
    ap.add_argument("--arc-mm", type=float, help="arc length for the bend/arc fallback")
    a = ap.parse_args(argv)
    table = export_map(a.results, a.design, a.material, a.radius_mm, a.arc_mm)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(table, indent=2) + "\n")
    print(f"wrote {a.out}: {len(table['pressure_kPa'])} levels")


if __name__ == "__main__":
    main()
