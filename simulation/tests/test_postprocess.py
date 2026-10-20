"""Post-solve path on fake FEBio logs: logs -> results.json, plots, renders.

The logs mimic FEBio's node_data file format (FECore/DataRecord.cpp, v4.13)
with a synthetic constant-curvature bend, so the whole read/metrics/render
path is exercised without FEBio installed.
"""

import json

import numpy as np
import pytest

from simulation.stage1_finger.analysis.metrics import synthetic_bend
from simulation.stage1_finger.backends import get_backend
from simulation.stage1_finger.backends.febio import LOG_FILES
from simulation.stage1_finger.mesh import fem_tags
from simulation.stage1_finger.run import postprocess
from simulation.stage1_finger.study import load_study

from conftest import STUDIES


def _write_log(path, name, times, rows_per_step):
    with open(path, "w") as f:
        for k, (t, rows) in enumerate(zip(times, rows_per_step), start=1):
            f.write(f"*Step  = {k}\n*Time  = {t:.9g}\n*Data  = {name}\n")
            for i, r in enumerate(rows, start=1):
                f.write(",".join([str(i)] + [f"{v:.12g}" for v in r]) + "\n")


@pytest.mark.parametrize("case", ["free", "blocked_force"])
def test_postprocess_fake_logs(smoke_mesh, tmp_path, case):
    mesh, _, _ = smoke_mesh
    study = load_study(STUDIES / "poc_smoke.toml")
    study.load.case = case
    (tmp_path / "mesh").mkdir()
    mesh.save(tmp_path / "mesh" / "mesh.npz")
    times = [0.25, 0.5, 0.75, 1.0]
    u = [synthetic_bend(mesh.nodes, 80 * t) for t in times]
    tip = mesh.node_sets["tip"]
    clamp = mesh.node_sets["clamp"]
    react = np.zeros((len(tip), 3))
    react[:, 2] = 0.01            # constraint pushes the finger toward +z
    _write_log(tmp_path / LOG_FILES["all_nodes"], "all_nodes", times, u)
    _write_log(tmp_path / LOG_FILES["tip"], "tip", times,
               [np.hstack([mesh.nodes[tip] + ui[tip], ui[tip], t * react]) for t, ui in zip(times, u)])
    _write_log(tmp_path / LOG_FILES["clamp"], "clamp", times,
               [np.zeros((len(clamp), 3)) for _ in times])

    rows = postprocess(tmp_path, mesh, study, get_backend("febio"), fem_tags(study.design),
                       "normal", render=(case == "free"))
    assert [r["pressure_kPa"] for r in rows] == pytest.approx([0, 0.5, 1.0, 1.5, 2.0])
    assert rows[-1]["bend_angle_deg"] == pytest.approx(80.0, abs=1e-6)
    data = json.loads((tmp_path / "results.json").read_text())
    assert data["reached_full_load"] and len(data["steps"]) == 5
    assert (tmp_path / "plots" / "p_vs_angle.png").exists()
    if case == "blocked_force":
        assert rows[-1]["blocked_force_N"] == pytest.approx(0.01 * len(tip))
    else:
        assert (tmp_path / "renders" / "deformed_100pct.png").exists()
