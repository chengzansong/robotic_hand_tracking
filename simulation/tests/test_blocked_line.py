"""blocked_line load case: a line support at the tip instead of the whole tip face.

The writer tests need no solve. The FEBio test solves a small tet10 cantilever
(48 elements, 12 x 4 x 3 mm, clamped at x = 0, pressure on the top face pushing
toward -z) twice, with the whole tip face held in z (blocked_force) and with
only the bottom edge of the tip face held (blocked_line). It checks that:
- unsupported tip nodes log about zero reaction, so summing the tip log
  (results.py) gives the support force for either case. FEBio logs a free
  node's internal-force sum, so a node that also carries pressure would log
  its share of the pressure load; the pressure therefore stops one cell short
  of the tip, as the finger's cavity does;
- the line carries z only: Rx = Ry = 0, and the reactions have no couple about
  the line (they act at one x within the line's spread);
- the tip rotates more on the line than on the face hold, which resists
  rotation once the face tilts.
Skipped unless febio_exe is configured.
"""

import itertools
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from simulation.core.materials import load_material
from simulation.stage1_finger.analysis.sweep import sweep_rows
from simulation.stage1_finger.backends.febio import LOG_FILES, febio_exe, run_febio, write_feb
from simulation.stage1_finger.backends.febio.results import log_status, read_node_data
from simulation.stage1_finger.backends.febio.writer import tip_line_nodes
from simulation.stage1_finger.mesh.fe_mesh import FEMesh, TET10_FACES, tet10_volumes
from simulation.stage1_finger.study import LoadCase, PartMaterial, load_study

from conftest import STUDIES

C1_STUDIES = {"cloth/c1_blocked_line_eco30_meshA.toml": "compare/sw6_blocked_eco30_meshA.toml",
              "cloth/c1_blocked_line_cloth_lo_meshA.toml": "cloth/cloth_lo_blocked_meshA.toml"}


def _write(mesh, tmp_path, case, materials=None, load=None, solver=None):
    path = tmp_path / f"{case}.feb"
    materials = materials or load_study(STUDIES / "poc_smoke.toml").materials
    info = write_feb(path, mesh, materials, load or LoadCase(case, 5.0, 4), bulk_to_shear=100.0,
                     contact={"enabled": False}, solver=solver)
    return ET.parse(path).getroot(), info


def _ids(root, name):
    el = root.find(f"Mesh/NodeSet[@name='{name}']")
    return None if el is None else [int(v) - 1 for v in el.text.replace("\n", "").split(",")]


def test_tip_line_node_set(half_mesh, tmp_path):
    mesh, _ = half_mesh
    root, info = _write(mesh, tmp_path, "blocked_line")
    line = _ids(root, "tip_line")
    tip = mesh.node_sets["tip"]
    assert sorted(line) == sorted(tip_line_nodes(mesh)) and info["tip_line_nodes"] == len(line)
    assert set(line) < set(tip) and len(line) >= 3
    p = mesh.nodes[line]
    assert np.allclose(p[:, 2], mesh.nodes[tip, 2].min())            # underside corner line
    assert np.allclose(p[:, 0], mesh.nodes[tip, 0].max())            # on the tip plane
    assert np.ptp(p[:, 1]) == pytest.approx(np.ptp(mesh.nodes[tip, 1]))   # across the width
    bcs = {b.get("name"): b for b in root.findall("Boundary/bc")}
    assert "tip_block" not in bcs and "symmetry" in bcs
    hold = bcs["tip_line_block"]
    assert hold.get("node_set") == "tip_line" and hold.get("type") == "prescribed displacement"
    assert hold.findtext("dof") == "z" and float(hold.findtext("value")) == 0.0
    # the tip log still covers the whole face, so results.py's summed reaction is unchanged
    assert root.find("Output/logfile/node_data[@name='tip']").get("node_set") == "tip"


@pytest.mark.parametrize("case", ["free", "blocked_force"])
def test_other_cases_write_no_tip_line(half_mesh, tmp_path, case):
    mesh, _ = half_mesh
    root, info = _write(mesh, tmp_path, case)
    assert _ids(root, "tip_line") is None and "tip_line_nodes" not in info
    assert root.find("Boundary/bc[@name='tip_line_block']") is None


@pytest.mark.parametrize("c1, source", C1_STUDIES.items())
def test_c1_studies_change_only_the_support(c1, source):
    a, b = load_study(STUDIES / c1), load_study(STUDIES / source)
    assert a.load.case == "blocked_line" and a.load.blocked and b.load.case == "blocked_force"
    ra, rb = dict(a.raw), dict(b.raw)
    assert ra.pop("name") == a.name == c1.split("/")[1][:-5]
    rb.pop("name")
    ra["load"] = {k: v for k, v in ra["load"].items() if k != "case"}
    rb["load"] = {k: v for k, v in rb["load"].items() if k != "case"}
    assert ra == rb


def test_sweep_reports_line_blocked_force(half_mesh):
    mesh, _ = half_mesh
    raw = [{"time": 1.0, "tip_reaction_N": [0.0, 0.0, -0.004], "clamp_reaction_N": [0.0] * 3}]
    rows, _, _ = sweep_rows(mesh, [1.0], np.zeros((1,) + mesh.nodes.shape), raw,
                            LoadCase("blocked_line", 5.0, 1))
    assert rows[-1]["blocked_force_N"] == pytest.approx(0.008)      # half model, doubled


# --- small FEBio solve ------------------------------------------------------

L, W, H = 12.0, 4.0, 3.0


def _beam(nx=4, ny=1, nz=2):
    """Structured tet10 box [0, L] x [0, W] x [0, H]: 6 Kuhn tets per cell."""
    g = np.stack(np.meshgrid(np.linspace(0, L, nx + 1), np.linspace(0, W, ny + 1),
                             np.linspace(0, H, nz + 1), indexing="ij"), -1)
    nodes = list(g.reshape(-1, 3))
    vid = np.arange(len(nodes)).reshape(nx + 1, ny + 1, nz + 1)
    tets = []
    for i, j, k in itertools.product(range(nx), range(ny), range(nz)):
        for perm in itertools.permutations(range(3)):
            path, c = [np.zeros(3, int)], np.zeros(3, int)
            for a in perm:
                c = c.copy()
                c[a] = 1
                path.append(c)
            tets.append([vid[i + p[0], j + p[1], k + p[2]] for p in path])
    mids, out = {}, []
    for t in tets:
        p = np.array([nodes[n] for n in t])
        if np.linalg.det(p[1:] - p[0]) < 0:
            t = [t[0], t[2], t[1], t[3]]
        e10 = list(t)
        for a, b in ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)):
            key = frozenset((t[a], t[b]))
            if key not in mids:
                mids[key] = len(nodes)
                nodes.append((nodes[t[a]] + nodes[t[b]]) / 2)
            e10.append(mids[key])
        out.append(e10)
    nodes, tets = np.array(nodes), np.array(out)
    faces = tets[:, TET10_FACES].reshape(-1, 6)
    keys = [frozenset(f[:3]) for f in faces]
    count = {k: keys.count(k) for k in set(keys)}
    outer = faces[[count[k] == 1 for k in keys]]
    on = lambda ax, v: outer[np.all(np.isclose(nodes[outer[:, :3], ax], v), axis=1)]
    top = on(2, H)
    top = top[nodes[top[:, :3], 0].max(1) <= L * (nx - 1) / nx + 1e-9]   # keep off the tip face
    surfaces = {"cavity": top, "tip": on(0, L)}
    node_sets = {"clamp": np.flatnonzero(np.isclose(nodes[:, 0], 0)),
                 "tip": np.unique(surfaces["tip"])}
    return FEMesh(nodes=nodes, parts={"top": tets}, surfaces=surfaces, node_sets=node_sets,
                  info={"fem_tags": {"bend_toward": "-z", "plane_tol_mm": 0.01}})


def test_beam_mesh_is_valid():
    mesh = _beam()
    assert mesh.n_elements == 48
    assert tet10_volumes(mesh.nodes, mesh.parts["top"]).sum() == pytest.approx(L * W * H)
    assert len(tip_line_nodes(mesh)) == 3            # two corners and the mid-node of one edge


@pytest.mark.skipif(febio_exe() is None, reason="febio_exe not configured")
def test_line_support_carries_no_moment(tmp_path):
    mesh = _beam()
    mat = {"top": PartMaterial(load_material("_placeholder_neo_hookean"), 0)}
    tip = mesh.node_sets["tip"]
    line = np.isin(tip, tip_line_nodes(mesh))
    res = {}
    for case in ("blocked_force", "blocked_line"):
        d = tmp_path / case
        d.mkdir()
        path = d / "model.feb"
        write_feb(path, mesh, mat, LoadCase(case, 1.0, 2), bulk_to_shear=100.0,
                  contact={"enabled": False}, solver={"threads": 1})
        info = run_febio(path, timeout_s=60, threads=1)
        assert info["returncode"] == 0 and log_status(d / "model.log") == "normal", case
        last = read_node_data(d / LOG_FILES["tip"])[-1]
        assert last["time"] == pytest.approx(1.0)
        res[case] = last["values"]                   # x;y;z;ux;uy;uz;Rx;Ry;Rz per tip node

    def couple(vals, held):     # moment about y of the held reactions, about their mean x
        x, rz = vals[held, 0], vals[held, 8]
        return np.sum((x - x.mean()) * rz)

    def tilt(vals):             # tip-face rotation: slope of deformed x over z
        return np.polyfit(vals[:, 2], vals[:, 0], 1)[0]

    v = res["blocked_line"]
    R = v[:, 6:9]
    F = -R[:, 2].sum()                               # finger pushes the support toward -z
    assert F > 0
    assert np.abs(R[~line]).max() < 1e-6 * F        # only the line is supported
    assert np.abs(R[:, :2]).max() < 1e-6 * F        # z only: free to slide in x
    # The face hold applies a couple once the face tilts (6e-3 F H here); the
    # line only up to its tiny spread in x across the width (1e-4 F H).
    face = res["blocked_force"]
    assert abs(couple(v, line)) < 1e-3 * F * H
    assert abs(couple(v, line)) < 0.1 * abs(couple(face, np.ones(len(tip), bool)))
    assert abs(tilt(v)) > abs(tilt(face)) > 0
