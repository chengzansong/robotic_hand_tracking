"""Run a Stage 1 study: mesh -> solver input -> solve -> metrics, plots, renders.

    python -m simulation.stage1_finger.run simulation/stage1_finger/studies/poc_smoke.toml
    python -m simulation.stage1_finger.run <study.toml> --no-solve    # mesh + .feb only
    python -m simulation.stage1_finger.run <study.toml> --estimate    # time estimate only

Before every solve the run probes the first pressure level(s) at full
settings (in <run>/estimate/), extrapolates the total wall time as a range,
and refuses to start the full solve if the upper bound exceeds --max-minutes
(default 20). A gated solve is killed after 1.5 x --max-minutes. Estimate and actual times go to metadata.json and to
simulation/outputs/estimates.csv, so the upper factor can be calibrated.

Writes simulation/outputs/<study>/<UTC timestamp>/:
    study.toml        copy of the study
    mesh/             mesh.msh (gmsh), mesh.vtu, mesh.npz (tags and node sets)
    model.feb         solver input
    model.log, *.txt  solver log and node_data logs (after a solve)
    results.json      metrics per pressure step (after a solve)
    displacements.npz nodal displacement per step, for re-rendering
    metadata.json     study, design, material fits used, mesh stats, git commit
    plots/, renders/  figures; renders can be regenerated with
                      python -m simulation.stage1_finger.render <run_dir>
"""

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from .backends import get_backend
from .mesh import build_mesh, fem_tags
from .study import load_study

OUTPUTS_DIR = Path(__file__).resolve().parents[1] / "outputs"
REPO_ROOT = Path(__file__).resolve().parents[2]
# Upper factor on the probe extrapolation. Calibrated from estimates.csv on
# 2026-10-20: converging runs took 0.90-1.20x the low estimate (5 runs); a run
# that struggled to a convergence failure took 2.4x, which the hard timeout
# (1.5 x --max-minutes) caps. Pass --upper-factor higher when the probe level is
# far less nonlinear than the rest of the sweep.
UPPER_FACTOR = 1.5


def git_commit():
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True,
                             text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO_ROOT,
                               capture_output=True, text=True).stdout.strip() != ""
        return {"commit": sha, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": "unknown", "dirty": None}


def _json(path, data):
    def default(o):
        if isinstance(o, (np.integer, np.floating)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, Path):
            return str(o)
        raise TypeError(type(o))
    Path(path).write_text(json.dumps(data, indent=2, default=default), encoding="utf-8")


def make_run_dir(study, out_root=OUTPUTS_DIR):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_dir = Path(out_root) / study.name / stamp
    run_dir.mkdir(parents=True, exist_ok=False)
    shutil.copy(study.path, run_dir / "study.toml")
    return run_dir


def estimate(run_dir, mesh, study, backend, tags, probe_levels=1, upper_factor=UPPER_FACTOR):
    """Solve the first `probe_levels` pressure levels and extrapolate the full run."""
    levels = study.load.levels()
    k = max(1, min(probe_levels, len(levels)))
    est_dir = Path(run_dir) / "estimate"
    est_dir.mkdir(exist_ok=True)
    written = backend.write(est_dir, mesh, study, bend_dof=tags["bend_toward"].lstrip("+-"),
                            end_time=study.load.fractions()[k - 1])
    info = backend.solve(est_dir, threads=int(written["solver"].get("threads", 0) or 0))
    probe_s = info["wall_time_s"]
    per_level = probe_s / k
    low = per_level * len(levels)
    out = {
        "probe_levels": k, "probe_kPa": levels[k - 1], "probe_s": probe_s,
        "probe_status": backend.status(est_dir), "probe_stats": info.get("stats", {}),
        "n_levels": len(levels), "upper_factor": upper_factor,
        "estimate_low_s": round(low, 1), "estimate_high_s": round(low * upper_factor, 1),
    }
    print(f"estimate: probe {k} level(s) to {levels[k - 1]:g} kPa took {probe_s:.0f} s "
          f"({out['probe_status']}); full run of {len(levels)} levels "
          f"~{low / 60:.1f} to {low * upper_factor / 60:.1f} min")
    return out


def _log_estimate(run_dir, study, est, actual_s=None, status=""):
    path = OUTPUTS_DIR / "estimates.csv"
    new = not path.exists()
    with open(path, "a", encoding="utf-8") as f:
        if new:
            f.write("run_dir,study,n_levels,probe_levels,probe_s,estimate_low_s,"
                    "estimate_high_s,actual_s,actual_over_low,status\n")
        ratio = round(actual_s / est["estimate_low_s"], 3) if actual_s and est["estimate_low_s"] else ""
        f.write(f"{run_dir},{study.name},{est['n_levels']},{est['probe_levels']},{est['probe_s']},"
                f"{est['estimate_low_s']},{est['estimate_high_s']},{actual_s or ''},{ratio},"
                f"{status}\n")


def run(study_path, solve=True, render=True, out_root=OUTPUTS_DIR, timeout_s=None,
        estimate_only=False, max_minutes=20.0, probe_levels=1, skip_estimate=False,
        upper_factor=UPPER_FACTOR):
    study = load_study(study_path)
    backend = get_backend(study.backend)
    run_dir = make_run_dir(study, out_root)
    print(f"run folder: {run_dir}")

    mesh, mesh_report = build_mesh(study.design, study.mesh, run_dir / "mesh")
    tags = fem_tags(study.design)
    s = mesh_report["stats"]
    print(f"mesh: {s['n_nodes']} nodes, {s['n_elements']} tet10, min SICN "
          f"{s['quality']['minSICN_min']}, cavity closed: {mesh_report['geometry']['cavity_closed']}")

    written = backend.write(run_dir, mesh, study, bend_dof=tags["bend_toward"].lstrip("+-"))
    meta = {
        "study": study.name,
        "study_file": str(study.path),
        "study_config": study.raw,
        "created_utc": run_dir.name,
        "design": study.design.id,
        "backend": study.backend,
        "materials": {part: {**pm.describe(), "solver": written["materials"][part]}
                      for part, pm in study.materials.items()},
        "bulk_to_shear": study.bulk_to_shear,
        "load": study.load.__dict__,
        "fem_tags": tags,
        "mesh": mesh_report,
        "solver_options": written["solver"],
        "contact": {"options": written["contact"], "slots": written["contact_slots"]},
        "levels_kPa": written["levels_kPa"],
        "git": git_commit(),
        "solve": {"status": "not run"},
    }
    _json(run_dir / "metadata.json", meta)
    if render:
        _render(run_dir)

    if not solve:
        print(f"--no-solve: wrote {written['path']}")
        return run_dir, 0
    if not backend.available():
        from .backends.febio.runner import INSTALL_HINT
        meta["solve"] = {"status": "skipped, solver not configured"}
        _json(run_dir / "metadata.json", meta)
        print(f"Solver not run. {INSTALL_HINT}\nRe-run without --no-solve once it is set up.",
              file=sys.stderr)
        return run_dir, 2

    est = None
    if not skip_estimate:
        est = estimate(run_dir, mesh, study, backend, tags, probe_levels=probe_levels,
                       upper_factor=upper_factor)
        meta["estimate"] = est
        _json(run_dir / "metadata.json", meta)
        if estimate_only:
            _log_estimate(run_dir, study, est, status="estimate only")
            return run_dir, 0
        if max_minutes and est["estimate_high_s"] > 60 * max_minutes:
            meta["solve"] = {"status": f"not run: estimate up to {est['estimate_high_s'] / 60:.1f} "
                                       f"min exceeds {max_minutes:g} min"}
            _json(run_dir / "metadata.json", meta)
            _log_estimate(run_dir, study, est, status="refused")
            print(meta["solve"]["status"], file=sys.stderr)
            return run_dir, 3
        if timeout_s is None and max_minutes:
            timeout_s = 1.5 * 60 * max_minutes

    threads = int(written["solver"].get("threads", 0) or 0)
    info = backend.solve(run_dir, timeout_s=timeout_s, threads=threads)
    info["status"] = backend.status(run_dir) if info["returncode"] != "timeout" else "timeout"
    meta["solve"] = info
    if est:
        est["actual_s"] = info["wall_time_s"]
        est["actual_over_low"] = round(info["wall_time_s"] / max(est["estimate_low_s"], 1e-9), 3)
        _log_estimate(run_dir, study, est, info["wall_time_s"], info["status"])
    _json(run_dir / "metadata.json", meta)
    print(f"solver finished: {info['status']} (return code {info['returncode']}, "
          f"{info['wall_time_s']} s)")

    rows = postprocess(run_dir, mesh, study, backend, tags, info["status"], render=render)
    if rows is None:
        return run_dir, 1
    ok = info["status"] == "normal" and rows[-1]["time"] > 1 - 1e-6
    return run_dir, 0 if ok else 1


def postprocess(run_dir, mesh, study, backend, tags, solver_status, render=True):
    """Solver logs -> results.json, displacements.npz, plots, renders. Returns rows."""
    from .analysis.curvature import turn_radius_sweep
    from .analysis.plots import plot_sweep, plot_turn_radius
    from .analysis.sweep import summarize, sweep_rows
    from .backends.febio.writer import CONTACT_DEFAULTS
    try:
        times, u, raw = backend.read(run_dir, len(mesh.nodes))
    except (OSError, ValueError) as err:
        print(f"could not read results: {err}", file=sys.stderr)
        return None
    tol = study.contact.get("penetration_tol_mm", CONTACT_DEFAULTS["penetration_tol_mm"])
    rows, u_all, p_all = sweep_rows(mesh, times, u, raw, study.load, axis=tags["axis"],
                                    bend_toward=tags["bend_toward"], penetration_tol_mm=tol)
    summary = summarize(rows, study.load, solver_status)
    if summary["penetration_warning"]:
        print(f"WARNING: slot walls interpenetrate by up to {-summary['min_slot_gap_mm']:.3f} mm "
              f"(from {summary['first_penetration_kPa']:g} kPa)", file=sys.stderr)
    if summary["curl_warning"]:
        print(f"WARNING: finger curls onto itself from {summary['first_curl_kPa']:g} kPa; "
              f"that contact is not modelled", file=sys.stderr)
    if not summary["reached_full_load"]:
        print(f"WARNING: stopped at {summary['last_converged_pressure_kPa']:g} of "
              f"{summary['requested_max_kPa']:g} kPa", file=sys.stderr)
    np.savez_compressed(run_dir / "displacements.npz", pressure_kPa=p_all, u=u_all,
                        time=np.array([r["time"] for r in rows]))
    turn = []
    if study.load.case == "free":
        try:
            turn, _ = turn_radius_sweep(mesh, u_all, p_all, axis=tags["axis"],
                                        bend_toward=tags["bend_toward"])
        except ValueError as err:
            print(f"turn radius skipped: {err}", file=sys.stderr)
    _json(run_dir / "results.json", {
        "study": study.name, "design": study.design.id, "load_case": study.load.case,
        **summary,
        "definitions": "see simulation/stage1_finger/analysis/metrics.py, contact.py and curvature.py",
        "steps": rows,
        "turn_radius": turn,
    })
    plot_sweep(rows, run_dir / "plots", title=f"{study.name} ({study.load.case})")
    if turn:
        plot_turn_radius(mesh, u_all, p_all, turn, run_dir / "plots", title=study.name,
                         axis=tags["axis"], bend_toward=tags["bend_toward"])
    if render:
        _render(run_dir)
    return rows


def _render(run_dir):
    from .analysis.render import render_run
    try:
        for path in render_run(run_dir):
            print(f"render: {path}")
    except Exception as err:     # a render failure must not lose a solve
        print(f"render failed: {err}", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Run a Stage 1 finger study.")
    ap.add_argument("study", help="path to a study TOML")
    ap.add_argument("--no-solve", action="store_true", help="mesh and write solver input only")
    ap.add_argument("--no-render", action="store_true", help="skip PyVista renders")
    ap.add_argument("--out-root", default=str(OUTPUTS_DIR))
    ap.add_argument("--timeout-s", type=float, default=None)
    ap.add_argument("--estimate", action="store_true", help="probe and print a time estimate only")
    ap.add_argument("--max-minutes", type=float, default=20.0,
                    help="refuse the full solve if the estimate's upper bound exceeds this")
    ap.add_argument("--probe-levels", type=int, default=1)
    ap.add_argument("--no-estimate", action="store_true", help="skip the estimate probe")
    ap.add_argument("--upper-factor", type=float, default=UPPER_FACTOR)
    args = ap.parse_args(argv)
    _, code = run(args.study, solve=not args.no_solve, render=not args.no_render,
                  out_root=args.out_root, timeout_s=args.timeout_s, estimate_only=args.estimate,
                  max_minutes=args.max_minutes, probe_levels=args.probe_levels,
                  skip_estimate=args.no_estimate, upper_factor=args.upper_factor)
    sys.exit(code)


if __name__ == "__main__":
    main()
