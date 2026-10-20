"""The .feb writer produces well-formed FEBio 4.0 XML with the right keywords."""

import xml.etree.ElementTree as ET

import pytest

from simulation.core.materials import list_materials
from simulation.stage1_finger.backends.febio import write_feb
from simulation.stage1_finger.study import LoadCase, PartMaterial, StudyError, load_study

from conftest import STUDIES

SECTIONS = ["Module", "Control", "Material", "Mesh", "MeshDomains", "Boundary", "Loads",
            "LoadData", "Output"]


def _write(mesh, tmp_path, part_materials, case="free", contact=None, load=None, solver=None):
    path = tmp_path / f"{case}.feb"
    info = write_feb(path, mesh, part_materials, load or LoadCase(case, 5.0, 4),
                     bulk_to_shear=100.0, contact=contact or {"enabled": False}, solver=solver)
    return ET.parse(path).getroot(), info, path


def test_smoke_feb_structure(smoke_mesh, tmp_path):
    mesh, _, _ = smoke_mesh
    study = load_study(STUDIES / "poc_smoke.toml")
    root, info, path = _write(mesh, tmp_path, study.materials)
    assert path.read_text(encoding="latin-1").startswith('<?xml version="1.0"')
    assert root.tag == "febio_spec" and root.get("version") == "4.0"
    assert [c.tag for c in root if isinstance(c.tag, str)] == SECTIONS
    assert root.find("Module").get("type") == "solid"
    assert root.findtext("Control/analysis") == "STATIC"

    mats = root.findall("Material/material")
    assert [m.get("type") for m in mats] == ["Mooney-Rivlin", "Mooney-Rivlin"]
    assert float(mats[0].findtext("c1")) == pytest.approx(0.040 / 2)   # mu/2 in MPa
    assert float(mats[0].findtext("c2")) == 0.0

    m = root.find("Mesh")
    assert len(m.findall("Nodes/node")) == len(mesh.nodes)
    elems = m.findall("Elements")
    assert [e.get("type") for e in elems] == ["tet10", "tet10"]
    assert sum(len(e) for e in elems) == mesh.n_elements
    assert all(len(e.text.split(",")) == 10 for e in elems[0][:20])
    tri6 = m.findall("Surface[@name='cavity']/tri6")
    assert len(tri6) == len(mesh.surfaces["cavity"])
    clamp = [int(v) for v in m.find("NodeSet[@name='clamp']").text.replace("\n", "").split(",")]
    assert sorted(clamp) == sorted(int(i) + 1 for i in mesh.node_sets["clamp"])

    doms = root.findall("MeshDomains/SolidDomain")
    assert {(d.get("name"), d.get("mat")) for d in doms} == {
        ("top", "top_material"), ("base", "base_material")}
    assert root.find("Loads/surface_load").get("type") == "pressure"
    assert float(root.findtext("Loads/surface_load/pressure")) == pytest.approx(0.005)
    bcs = root.findall("Boundary/bc")
    # prescribed zero displacements, so FEBio logs clamp reactions
    assert [(b.get("name"), b.get("type"), b.findtext("dof")) for b in bcs] == [
        (f"clamp_{d}", "prescribed displacement", d) for d in "xyz"]
    assert len(root.findall("LoadData/load_controller[@id='2']/points/pt")) == 4
    assert {n.get("name") for n in root.findall("Output/logfile/node_data")} == {
        "tip", "clamp", "all_nodes"}


def test_blocked_force_bc(smoke_mesh, tmp_path):
    mesh, _, _ = smoke_mesh
    study = load_study(STUDIES / "poc_smoke.toml")
    root, _, _ = _write(mesh, tmp_path, study.materials, case="blocked_force")
    tip = root.find("Boundary/bc[@name='tip_block']")
    assert tip.get("node_set") == "tip"
    assert tip.get("type") == "prescribed displacement"
    assert tip.findtext("dof") == "z"
    assert float(tip.findtext("value")) == 0.0


@pytest.mark.skipif("ecoflex_00-30" not in list_materials(), reason="card not present")
def test_real_card_keywords(smoke_mesh, tmp_path):
    from simulation.core.materials import load_material
    mesh, _, _ = smoke_mesh
    eco = load_material("ecoflex_00-30")
    by_model = {f.model: i for i, f in enumerate(eco.fits)}
    expect = {"yeoh": "Yeoh", "ogden": "Ogden", "neo_hookean": "Mooney-Rivlin"}
    for model, idx in by_model.items():
        pm = {"top": PartMaterial(eco, idx), "base": PartMaterial(eco, 0)}
        root, _, _ = _write(mesh, tmp_path, pm)
        assert root.find("Material/material[@id='1']").get("type") == expect[model]


def test_study_fit_selection_and_missing_card(tmp_path):
    s = tmp_path / "s.toml"
    s.write_text('design = "finger_poc"\n[materials]\n'
                 'top = { material = "_placeholder_neo_hookean", fit = 0 }\n'
                 'base = "_placeholder_neo_hookean"\n[load]\npressure_kPa = 1.0\nsteps = 2\n')
    study = load_study(s)
    assert study.materials["top"].fit_index == 0
    assert study.materials["top"].describe()["model"] == "neo_hookean"
    s.write_text('design = "finger_poc"\n[materials]\ntop = "no_such_card"\n')
    with pytest.raises(StudyError, match="no_such_card.*not found"):
        load_study(s)
    s.write_text('design = "finger_poc"\n[materials]\n'
                 'top = { material = "_placeholder_neo_hookean", fit = 5 }\n')
    with pytest.raises(StudyError, match="out of range"):
        load_study(s)


def test_contact_section(smoke_mesh, tmp_path):
    mesh, _, _ = smoke_mesh
    study = load_study(STUDIES / "poc_smoke.toml")
    root, info, _ = _write(mesh, tmp_path, study.materials, contact={"enabled": True})
    tags = [c.tag for c in root if isinstance(c.tag, str)]
    assert tags.index("Loads") < tags.index("Contact") < tags.index("LoadData")
    pairs = root.findall("Mesh/SurfacePair")
    assert len(pairs) == 10 == len(info["contact_slots"])
    surfaces = {sf.get("name") for sf in root.findall("Mesh/Surface")}
    for pair in pairs:
        assert {pair.findtext("primary"), pair.findtext("secondary")} <= surfaces
    contacts = root.findall("Contact/contact")
    assert {c.get("type") for c in contacts} == {"sliding-elastic"}
    assert {c.get("surface_pair") for c in contacts} == {p.get("name") for p in pairs}
    c0 = contacts[0]
    assert c0.findtext("laugon") == "PENALTY" and c0.findtext("two_pass") == "1"
    assert c0.findtext("auto_penalty") == "1" and float(c0.findtext("fric_coeff")) == 0.0
    assert "contact pressure" in {v.get("type") for v in root.findall("Output/plotfile/var")}


def test_levels_and_solver_options(smoke_mesh, tmp_path):
    mesh, _, _ = smoke_mesh
    study = load_study(STUDIES / "poc_smoke.toml")
    load = LoadCase(levels_kPa=[0.5, 1.0, 2.0, 5.0])
    root, info, _ = _write(mesh, tmp_path, study.materials, load=load,
                           solver={"symmetric_stiffness": "non-symmetric", "max_ups": 0,
                                   "pressure_symmetric": False})
    pts = [tuple(map(float, p.text.split(",")))
           for p in root.findall("LoadData/load_controller[@id='2']/points/pt")]
    assert pts == pytest.approx([(0.1, 0.1), (0.2, 0.1), (0.4, 0.2), (1.0, 0.6)])
    assert float(root.findtext("Loads/surface_load/pressure")) == pytest.approx(0.005)
    assert root.findtext("Loads/surface_load/symmetric_stiffness") == "0"
    assert root.findtext("Control/solver/symmetric_stiffness") == "non-symmetric"
    assert root.findtext("Control/solver/qn_method/max_ups") == "0"
    n, dt = int(root.findtext("Control/time_steps")), float(root.findtext("Control/step_size"))
    assert n * dt == pytest.approx(1.0)


def test_half_model_symmetry_bc(half_mesh, tmp_path):
    mesh, _ = half_mesh
    study = load_study(STUDIES / "poc_smoke.toml")
    root, _, _ = _write(mesh, tmp_path, study.materials)
    bc = root.find("Boundary/bc[@name='symmetry']")
    assert bc.get("type") == "zero displacement"
    assert [bc.findtext(f"{d}_dof") for d in "xyz"] == ["0", "1", "0"]
    ids = root.find("Mesh/NodeSet[@name='symmetry']").text.replace("\n", "").split(",")
    assert len(ids) == len(mesh.node_sets["symmetry"])


def test_study_levels_validation(tmp_path):
    s = tmp_path / "s.toml"
    s.write_text('design = "finger_poc"\n[load]\nlevels_kPa = [0.5, 1, 3]\n')
    st = load_study(s)
    assert st.load.max_kPa == 3.0 and st.load.steps == 3
    s.write_text('design = "finger_poc"\n[load]\nlevels_kPa = [1, 0.5]\n')
    with pytest.raises(StudyError, match="increasing"):
        load_study(s)
