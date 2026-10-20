"""The time estimate and the 20-minute gate, with a fake backend (no FEBio)."""

import types

import pytest

from simulation.stage1_finger import run as run_mod
from simulation.stage1_finger.study import load_study

from conftest import STUDIES


def _fake_backend(probe_s, calls):
    def write(run_dir, mesh, study, bend_dof, end_time=1.0):
        calls.append(("write", run_dir.name, end_time))
        return {"path": "x", "materials": {p: {} for p in study.materials}, "solver": {},
                "contact": {}, "contact_slots": [], "levels_kPa": study.load.levels()}

    def solve(run_dir, timeout_s=None, threads=0):
        calls.append(("solve", run_dir.name))
        return {"wall_time_s": probe_s, "returncode": 0, "stats": {}}

    return types.SimpleNamespace(write=write, solve=solve, status=lambda d: "normal",
                                 available=lambda: True, read=None)


@pytest.mark.parametrize("probe_s, expect_run", [(10.0, True), (200.0, False)])
def test_estimate_gate(monkeypatch, tmp_path, smoke_mesh, probe_s, expect_run):
    mesh, report, _ = smoke_mesh
    calls = []
    monkeypatch.setattr(run_mod, "get_backend", lambda name: _fake_backend(probe_s, calls))
    monkeypatch.setattr(run_mod, "build_mesh", lambda design, spec, out: (mesh, report))
    monkeypatch.setattr(run_mod, "OUTPUTS_DIR", tmp_path)
    monkeypatch.setattr(run_mod, "postprocess", lambda *a, **k: [{"time": 1.0}])
    study = load_study(STUDIES / "poc_contact_check.toml")     # 6 levels
    run_dir, code = run_mod.run(study.path, render=False, out_root=tmp_path, max_minutes=20)
    est = __import__("json").loads((run_dir / "metadata.json").read_text())["estimate"]
    assert est["estimate_low_s"] == pytest.approx(6 * probe_s)
    assert est["estimate_high_s"] == pytest.approx(6 * probe_s * run_mod.UPPER_FACTOR)
    probe = [c for c in calls if c[0] == "write" and c[1] == "estimate"]
    assert probe and probe[0][2] == pytest.approx(0.5 / 3.0)    # first level of 3 kPa
    full_solves = [c for c in calls if c == ("solve", run_dir.name)]
    assert bool(full_solves) == expect_run
    assert code == (0 if expect_run else 3)
    assert (tmp_path / "estimates.csv").exists()
