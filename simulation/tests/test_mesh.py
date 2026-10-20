"""The finger_poc mesh: parts, tags, orientation, volumes, and extents."""

import numpy as np
import pytest

from simulation.stage1_finger.mesh.fe_mesh import (
    TET10_FACES, FEMesh, boundary_loops, cavity_volume, enclosed_volume, tet10_jacobians, tet10_volumes,
    tri6_to_linear)

EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]


def test_two_parts_and_tags(smoke_mesh):
    mesh, report, _ = smoke_mesh
    assert set(mesh.parts) == {"top", "base"}
    assert all(len(t) > 0 for t in mesh.parts.values())
    for name in ("cavity", "clamp", "tip"):
        assert len(mesh.surfaces[name]) > 0
        assert len(mesh.node_sets[name]) > 0
    assert len(mesh.surfaces["interface"]) > 0
    g = report["geometry"]
    assert g["cavity_closed"] and g["n_shells"] == 2
    assert g["overlap"]["volume_mm3"] == pytest.approx(338.0, abs=1e-3)


def test_extents_and_volumes(smoke_mesh):
    mesh, report, _ = smoke_mesh
    assert np.ptp(mesh.nodes, axis=0) == pytest.approx([112.0, 15.0, 19.0], abs=1e-4)
    tets = mesh.all_tets()
    vol = tet10_volumes(mesh.nodes, tets)
    assert vol.min() > 0
    assert (tet10_jacobians(mesh.nodes, tets) > 0).all()
    # top 12340 + base 6720 - 338 overlap
    assert vol.sum() == pytest.approx(18722.0, rel=1e-6)
    assert report["stats"]["quality"]["minSICN_min"] > 0.05


def test_tet10_midnodes_in_vtk_febio_order(smoke_mesh):
    mesh, _, _ = smoke_mesh
    t = mesh.all_tets()[:200]
    x = mesh.nodes
    for k, (i, j) in enumerate(EDGES, start=4):
        assert np.allclose(x[t[:, k]], (x[t[:, i]] + x[t[:, j]]) / 2, atol=1e-6)


def test_cavity_volume_and_orientation(smoke_mesh):
    mesh, _, _ = smoke_mesh
    v, n_open = cavity_volume(mesh.nodes, mesh.surfaces["cavity"])
    # 11 chambers 6 x 11 x 13 mm plus 10 feed slots 4 x 2 x 2 mm
    assert n_open == 0
    assert v == pytest.approx(11 * 6 * 11 * 13 + 10 * 16, rel=1e-6)
    # facets point out of the solid: their normals point away from the nearest tet
    tri = mesh.surfaces["cavity"][:50]
    p = mesh.nodes[tri[:, :3]]
    n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    c_tet = mesh.nodes[mesh.all_tets()[:, :4]].mean(1)
    for normal, centroid in zip(n, p.mean(1)):
        nearest = c_tet[np.argmin(np.linalg.norm(c_tet - centroid, axis=1))]
        assert normal @ (centroid - nearest) > 0


def test_clamp_and_tip_planes(smoke_mesh):
    mesh, _, _ = smoke_mesh
    assert np.allclose(mesh.nodes[mesh.node_sets["clamp"], 0], 0.0, atol=1e-6)
    assert np.allclose(mesh.nodes[mesh.node_sets["tip"], 0], 112.0, atol=1e-6)
    # each end face covers the full 15 x 19 mm section
    for name in ("clamp", "tip"):
        t = mesh.nodes[tri6_to_linear(mesh.surfaces[name])]
        area = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1).sum() / 2
        assert area == pytest.approx(15 * 19, rel=1e-6)


def test_save_load_roundtrip(smoke_mesh, tmp_path):
    mesh, _, _ = smoke_mesh
    mesh.save(tmp_path / "m.npz")
    back = FEMesh.load(tmp_path / "m.npz")
    assert np.array_equal(back.nodes, mesh.nodes)
    assert np.array_equal(back.surfaces["cavity"], mesh.surfaces["cavity"])
    assert back.info["geometry"]["cavity_closed"]


def test_open_surface_is_capped():
    # unit cube surface with outward normals, top face (z = 1) removed
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    quads = [(0, 3, 2, 1), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    tris = np.array([t for a, b, c, d in quads for t in ((a, b, c), (a, c, d))])
    loops = boundary_loops(tris)
    assert len(loops) == 1 and len(loops[0]) == 4
    assert enclosed_volume(pts, tris, caps=loops) == pytest.approx(1.0)


def test_tet10_face_table_points_outward():
    corners = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    x = np.vstack([corners] + [(corners[i] + corners[j]) / 2 for i, j in EDGES])
    for f in TET10_FACES:
        p = x[f[:3]]
        assert np.cross(p[1] - p[0], p[2] - p[0]) @ (p.mean(0) - corners.mean(0)) > 0
        for k, (a, b) in zip(f[3:], ((0, 1), (1, 2), (2, 0))):
            assert np.allclose(x[k], (p[a] + p[b]) / 2)


def test_slot_walls_found(smoke_mesh):
    mesh, report, _ = smoke_mesh
    slots = report["geometry"]["slots"]
    assert len(slots) == 10
    assert [s["lo_mm"] for s in slots] == pytest.approx([10.0 * k for k in range(1, 11)])
    assert all(s["gap_mm"] == pytest.approx(2.0) for s in slots)
    for s in mesh.info["slots"]:
        lo = mesh.nodes[np.unique(mesh.surfaces[f"{s}_lo"])]
        hi = mesh.nodes[np.unique(mesh.surfaces[f"{s}_hi"])]
        assert np.ptp(lo[:, 0]) < 1e-6 and np.ptp(hi[:, 0]) < 1e-6
        assert hi[0, 0] - lo[0, 0] == pytest.approx(2.0)


def test_slot_gaps_undeformed(smoke_mesh):
    from simulation.stage1_finger.analysis.contact import slot_gaps
    mesh, _, _ = smoke_mesh
    gaps = slot_gaps(mesh, mesh.nodes)
    assert len(gaps) == 10 and all(g == pytest.approx(2.0, abs=1e-6) for g in gaps.values())
    # push every hi wall node 2.3 mm toward -x: walls overlap by 0.3 mm
    x = mesh.nodes.copy()
    for s in mesh.info["slots"]:
        x[np.unique(mesh.surfaces[f"{s}_hi"]), 0] -= 2.3
    assert min(slot_gaps(mesh, x).values()) == pytest.approx(-0.3, abs=1e-6)


def test_half_model_matches_full(half_mesh, smoke_mesh):
    from simulation.stage1_finger.mesh.fe_mesh import mirror_half
    half, report = half_mesh
    full_ref, _, _ = smoke_mesh
    assert report["stats"]["solid_volume_mm3"] == pytest.approx(18722.0 / 2, rel=1e-6)
    assert report["geometry"]["cavity_volume_mm3"] == pytest.approx(9598.0, rel=1e-6)
    assert np.allclose(half.nodes[half.node_sets["symmetry"], 1], 7.5)
    full = mirror_half(half)
    v, _ = cavity_volume(full.nodes, full.surfaces["cavity"])
    assert v == pytest.approx(9598.0, rel=1e-6)
    assert tet10_volumes(full.nodes, full.all_tets()).min() > 0
    assert np.ptp(full.nodes, axis=0) == pytest.approx(np.ptp(full_ref.nodes, axis=0))


def test_wall_size_from_thinnest_wall():
    from simulation.core.designs import list_designs, load_design
    from simulation.stage1_finger.mesh.build import MeshSpec, resolve_wall_size, thinnest_wall_mm
    poc = load_design("finger_poc")
    assert thinnest_wall_mm(poc) == 1.0
    assert resolve_wall_size(poc, MeshSpec(wall_size_mm=1.5)) == 1.5
    assert resolve_wall_size(poc, MeshSpec(wall_size_mm=1.5, wall_elems=2)) == 0.5
    if "pn_base4_wall0p75" in list_designs():
        thin = load_design("pn_base4_wall0p75")
        assert thinnest_wall_mm(thin) == 0.75
        assert resolve_wall_size(thin, MeshSpec(wall_size_mm=1.0, wall_elems=1.5)) == 0.5
