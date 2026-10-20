"""Metrics, the log parser, and the render path, on synthetic data (no FEBio)."""

import numpy as np
import pytest

from simulation.stage1_finger.analysis.metrics import step_metrics, synthetic_bend
from simulation.stage1_finger.analysis.sweep import sweep_rows
from simulation.stage1_finger.backends.febio.results import read_node_data
from simulation.stage1_finger.mesh.fe_mesh import principal_stretches
from simulation.stage1_finger.study import LoadCase


@pytest.mark.parametrize("angle", [0.0, 30.0, 120.0])
def test_bend_angle_of_synthetic_bend(smoke_mesh, angle):
    mesh, _, _ = smoke_mesh
    m = step_metrics(mesh, synthetic_bend(mesh.nodes, angle))
    assert m["bend_angle_deg"] == pytest.approx(angle, abs=1e-6)
    assert abs(m["out_of_plane_deg"]) < 1e-6


def test_sweep_rows_unwrap_past_180(smoke_mesh):
    mesh, _, _ = smoke_mesh
    angles = [60, 120, 180, 240]
    u = np.stack([synthetic_bend(mesh.nodes, a) for a in angles])
    rows, _, p = sweep_rows(mesh, np.array([0.25, 0.5, 0.75, 1.0]), u,
                            [{"time": t} for t in (0.25, 0.5, 0.75, 1.0)],
                            LoadCase("free", 20.0, 4))
    assert [r["bend_angle_deg"] for r in rows] == pytest.approx([0] + angles, abs=1e-6)
    assert p.tolist() == [0, 5, 10, 15, 20]


def test_uniform_stretch():
    c = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float) * [2.0, 1.0, 3.0]
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    nodes = np.vstack([c] + [(c[i] + c[j]) / 2 for i, j in edges])
    tets = np.arange(10)[None, :]
    lam = principal_stretches(nodes, nodes * [2.0, 0.5, 1.0], tets)
    assert lam[0] == pytest.approx([0.5, 1.0, 2.0])


def test_node_data_parser(tmp_path):
    f = tmp_path / "tip.txt"
    f.write_text("*Step  = 1\n*Time  = 0.25\n*Data  = tip\n1,0,0,1\n2,0,1,2\n"
                 "*Step  = 2\n*Time  = 0.5\n*Data  = tip\n1,0,0,3\n2,0,1,4\n")
    recs = read_node_data(f)
    assert [r["time"] for r in recs] == [0.25, 0.5]
    assert recs[1]["values"].tolist() == [[0, 0, 3], [0, 1, 4]]


def test_render_run_on_synthetic_sweep(smoke_mesh, tmp_path):
    pytest.importorskip("pyvista")
    from simulation.stage1_finger.analysis.render import render_run
    mesh, _, _ = smoke_mesh
    (tmp_path / "mesh").mkdir()
    mesh.save(tmp_path / "mesh" / "mesh.npz")
    p = np.linspace(0, 20, 5)
    np.savez(tmp_path / "displacements.npz", pressure_kPa=p,
             u=np.stack([synthetic_bend(mesh.nodes, 4.5 * q) for q in p]))
    written = render_run(tmp_path)
    names = {w.name for w in written}
    assert {"mesh_tags.png", "deformed_100pct.png", "section_max.png",
            "sweep_montage.png"} <= names
    assert all(w.stat().st_size > 10_000 for w in written)


def test_frame_indices_distinct():
    from simulation.stage1_finger.analysis.render import _frame_indices
    assert _frame_indices([0, 1, 2, 4]) == [1, 2, 3]          # was [1, 2, 2, 3]
    assert _frame_indices([0, 0.5, 1, 2, 3, 4, 5]) == [2, 3, 5, 6]


def test_center_sections_on_synthetic_half_model(half_mesh, tmp_path):
    pytest.importorskip("pyvista")
    from PIL import Image
    from simulation.stage1_finger.analysis.render import render_center_sections
    mesh, _ = half_mesh
    (tmp_path / "mesh").mkdir()
    mesh.save(tmp_path / "mesh" / "mesh.npz")
    p = np.array([0.0, 1.0, 2.0])
    np.savez(tmp_path / "displacements.npz", pressure_kPa=p,
             u=np.stack([synthetic_bend(mesh.nodes, 30 * q) for q in p]))
    paths, summary = render_center_sections(tmp_path)
    out = tmp_path / "renders" / "center"
    pngs = sorted(out.glob("center_*.png"))
    assert len(pngs) == 2 * len(p)
    assert all(f.stat().st_size > 10_000 for f in pngs)
    for view in ("full", "zoom"):
        with Image.open(out / f"center_{view}.gif") as gif:
            assert gif.n_frames == len(p)
    assert [s["pressure_kPa"] for s in summary] == p.tolist()
    assert all(s["inverted"] == 0 for s in summary)
    assert set(paths) >= set(pngs)
