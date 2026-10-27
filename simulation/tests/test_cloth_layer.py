"""Cloth (strain-limiting fabric) layer: generator, mesh tags and bonding, study and .feb.

A cloth design is generated into a temporary models folder (generated designs
are gitignored), meshed as a coarse half model, and checked for: three parts
matched by volume; a cloth slab of the right volume that shares its nodes with
the base on its lower, upper, and side faces (bonded, conformal); the cloth end
faces in the clamp and tip sets; the same total solid and cavity as finger_poc.
"""

import tomllib
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from simulation.core.designs import load_design
from simulation.core.materials import list_materials, load_material
from simulation.models.tools import make_pneunet as mp
from simulation.stage1_finger.backends.febio import write_feb
from simulation.stage1_finger.mesh import MeshSpec, build_mesh
from simulation.stage1_finger.mesh.fe_mesh import tet10_volumes, tri6_to_linear
from simulation.stage1_finger.study import LoadCase, PartMaterial, StudyError, load_study

Z0, T = 2.0, 1.0      # layer z = 2..3 mm, inset 1 mm: 112 x 13 x 1 mm


@pytest.fixture(scope="module")
def cloth_design(tmp_path_factory):
    models = tmp_path_factory.mktemp("models")
    old = mp.MODELS_DIR
    mp.MODELS_DIR = models
    try:
        vols = mp.make("t_cloth", mp.PneuNet(cloth_t_mm=T, cloth_z_mm=Z0))
    finally:
        mp.MODELS_DIR = old
    return load_design("t_cloth", models_dir=models), vols, models


@pytest.fixture(scope="module")
def cloth_mesh(cloth_design):
    design, _, _ = cloth_design
    return build_mesh(design, MeshSpec(size_mm=4.0, wall_size_mm=1.5, symmetry="half_y"))


def test_generator_volumes_and_design(cloth_design):
    design, vols, models = cloth_design
    assert vols["cloth_mm3"] == pytest.approx(112 * 13 * T, rel=1e-9)
    assert vols["base_mm3"] == pytest.approx(112 * 15 * 4 - 112 * 13 * T, rel=1e-9)
    assert vols["top_mm3"] == pytest.approx(12340.0, rel=1e-6)
    assert [p.name for p in design.parts] == ["top", "base", "cloth"]
    cloth = design.part("cloth")
    assert cloth.layer_thickness_mm == T and (cloth.warp_axis, cloth.weft_axis) == ("x", "y")
    assert cloth.default_material == "fabric_woven_cotton"
    with open(models / "t_cloth" / "design.toml", "rb") as f:
        par = tomllib.load(f)["parametric"]
    assert (par["cloth_t_mm"], par["cloth_z_mm"], par["cloth_inset_mm"]) == (T, Z0, 1.0)


def test_generator_defaults_unchanged_and_checks():
    p = mp.PneuNet()
    assert not p.has_cloth
    # the cloth keys are left out of design.toml [parametric] when there is no cloth
    assert set(mp.CLOTH_KEYS) <= set(p.__dict__)
    with pytest.raises(ValueError, match="below the lip"):
        mp.PneuNet(cloth_t_mm=1.0, cloth_z_mm=2.5).check()      # 2.5..3.5 hits the pegs
    with pytest.raises(ValueError, match="inset"):
        mp.PneuNet(cloth_t_mm=1.0, cloth_z_mm=1.0, cloth_inset_mm=0.0).check()
    mp.PneuNet(cloth_t_mm=1.0, cloth_z_mm=2.0).check()          # top face at the lip, allowed


def test_three_parts_bonded(cloth_mesh):
    mesh, report = cloth_mesh
    assert set(mesh.parts) == {"top", "base", "cloth"}
    g = report["geometry"]
    assert g["n_shells"] == 2      # half model: the cavity is open on the symmetry plane
    assert g["overlap"]["volume_mm3"] == pytest.approx(338.0 / 2, abs=1e-3)   # lip and pegs only
    assert g["cavity_volume_mm3"] == pytest.approx(9598.0, rel=1e-6)
    vol = {k: tet10_volumes(mesh.nodes, v).sum() for k, v in mesh.parts.items()}
    assert vol["cloth"] == pytest.approx(112 * 6.5 * T, rel=1e-6)
    assert sum(vol.values()) == pytest.approx(18722.0 / 2, rel=1e-6)   # same solid as finger_poc

    cloth_nodes = np.unique(mesh.parts["cloth"])
    base_nodes = set(np.unique(mesh.parts["base"]).tolist())
    # the lip and pegs (top material) sit on part of the upper face
    solid_nodes = base_nodes | set(np.unique(mesh.parts["top"]).tolist())
    x = mesh.nodes[cloth_nodes]
    assert x[:, 2].min() == pytest.approx(Z0) and x[:, 2].max() == pytest.approx(Z0 + T)
    assert x[:, 1].min() == pytest.approx(1.0) and x[:, 1].max() == pytest.approx(7.5)
    # every cloth node on its lower and side faces is also a base node, and on its
    # upper face a base or lip/peg node: bonded
    for name, on_face in (("lower", np.isclose(x[:, 2], Z0)), ("upper", np.isclose(x[:, 2], Z0 + T)),
                          ("side", np.isclose(x[:, 1], 1.0))):
        ids = cloth_nodes[on_face]
        assert len(ids) > 20, name
        assert all(int(i) in (solid_nodes if name == "upper" else base_nodes) for i in ids), name
    # no cloth node strictly inside the layer is shared (no node merging across it)
    inner = cloth_nodes[(x[:, 2] > Z0 + 1e-6) & (x[:, 2] < Z0 + T - 1e-6) & (x[:, 1] > 1 + 1e-6)
                        & (x[:, 1] < 7.5 - 1e-6) & (x[:, 0] > 1e-6) & (x[:, 0] < 112 - 1e-6)]
    assert not any(int(i) in solid_nodes for i in inner)
    # the bonded faces are in the interface group: lower + upper + side faces of the half layer
    tri = mesh.nodes[tri6_to_linear(mesh.surfaces["interface"])]
    cz = tri[..., 2].mean(1)
    cy = tri[..., 1].mean(1)
    area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    on_cloth = ((np.isclose(cz, Z0) | np.isclose(cz, Z0 + T)) & (cy > 1 - 1e-6)) | (
        np.isclose(cy, 1.0) & (cz > Z0 - 1e-6) & (cz < Z0 + T + 1e-6))
    assert area[on_cloth].sum() == pytest.approx(2 * 112 * 6.5 + 112 * T, rel=1e-6)
    # the cloth's end faces are clamped and part of the tip
    for name, xv in (("clamp", 0.0), ("tip", 112.0)):
        s = set(mesh.node_sets[name].tolist())
        ends = cloth_nodes[np.isclose(x[:, 0], xv)]
        assert len(ends) > 3 and all(int(i) in s for i in ends)


def test_layer_size_refines_only_the_cloth(cloth_design, cloth_mesh):
    design, _, _ = cloth_design
    plain, _ = cloth_mesh
    fine, rep = build_mesh(design, MeshSpec(size_mm=4.0, wall_size_mm=1.5, symmetry="half_y",
                                            part_size_mm={"cloth": 1.0}))
    assert len(fine.parts["cloth"]) > 3 * len(plain.parts["cloth"])
    assert rep["stats"]["quality"]["minSICN_min"] > 0.05
    with pytest.raises(ValueError, match="part_size_mm"):
        build_mesh(design, MeshSpec(size_mm=4.0, wall_size_mm=1.5, part_size_mm={"nope": 1.0}))


def test_sheet_card_kind():
    assert "fabric_woven_cotton" not in list_materials()            # solids only by default
    assert "fabric_woven_cotton" in list_materials(kind="sheet")
    m = load_material("fabric_woven_cotton")
    assert m.kind == "sheet" and m.default_fit.is_sheet
    assert m.default_fit.warp_weft_N_per_mm() == (15.0, 10.0)
    with pytest.raises(ValueError):
        m.default_fit.initial_shear_modulus_kPa()


def _study(tmp_path, models, body):
    path = tmp_path / "s.toml"
    path.write_text('design = "t_cloth"\n' + body, encoding="utf-8")
    import simulation.stage1_finger.study as st
    old = st.load_design
    st.load_design = lambda d: load_design(d, models_dir=models)
    try:
        return load_study(path)
    finally:
        st.load_design = old


def test_study_and_feb(cloth_design, cloth_mesh, tmp_path):
    _, _, models = cloth_design
    mesh, _ = cloth_mesh
    s = _study(tmp_path, models, '[materials]\nbase = "ecoflex_00-30"\n'
               'cloth = { material = "fabric_woven_cotton", stiffness_scale = 2.0 }\n')
    pm = s.materials["cloth"]
    assert pm.is_sheet and pm.matrix.material.id == "ecoflex_00-30" and pm.stiffness_scale == 2.0
    assert pm.layer == {"thickness_mm": T, "warp_axis": "x", "weft_axis": "y"}
    assert pm.describe()["matrix"]["material"] == "ecoflex_00-30"
    path = tmp_path / "c.feb"
    info = write_feb(path, mesh, s.materials, LoadCase("free", 2.0, 2), 100.0,
                     contact={"enabled": False})
    root = ET.parse(path).getroot()
    mats = {m.get("name"): m for m in root.findall("Material/material")}
    assert mats["top_material"].get("type") == "Yeoh"
    cl = mats["cloth_material"]
    assert cl.get("type") == "solid mixture"
    sol = cl.findall("solid")
    assert [x.get("type") for x in sol] == ["neo-Hookean", "fiber-pow-linear", "fiber-pow-linear"]
    assert [x.findtext("fiber") for x in sol[1:]] == ["1,0,0", "0,1,0"]
    assert sol[1].find("fiber").get("type") == "vector"
    assert info["materials"]["cloth"]["params_MPa"]["warp_N_per_mm"] == pytest.approx(30.0)
    doms = {(d.get("name"), d.get("mat")) for d in root.findall("MeshDomains/SolidDomain")}
    assert ("cloth", "cloth_material") in doms

    # plain silicone on the cloth part (a geometry control) is a normal solid
    s2 = _study(tmp_path, models, '[materials]\ncloth = "ecoflex_00-30"\n')
    assert not s2.materials["cloth"].is_sheet
    with pytest.raises(StudyError, match="sheet materials only"):
        _study(tmp_path, models, '[materials]\ncloth = { material = "ecoflex_00-30", matrix = "ecoflex_00-30" }\n')
    with pytest.raises(StudyError, match="layer_thickness_mm"):
        _study(tmp_path, models, '[materials]\nbase = "fabric_woven_cotton"\n')
