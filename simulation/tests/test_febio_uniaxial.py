"""Single tet10 in uniaxial tension, solved by FEBio, against analytic neo-Hookean.

The element is the corner tet (0,0,0) (1,0,0) (0,1,0) (0,0,1) mm. Symmetry
planes hold ux on x = 0, uy on y = 0, and uz on z = 0. Every other node gets
ux = (lam - 1) x with uy and uz free, so the exact solution is homogeneous
uniaxial stretch, which tet10 represents exactly. The x = 0 face has reference
area 0.5 mm^2 and carries 0.5 * P, with nominal stress P = mu (lam - lam^-2)
for an incompressible neo-Hookean solid. Only 2/3 of that shows up as nodal
reaction on the x = 0 set: under uniform stress a tri6 face loads its three
mid-side nodes equally, and mid-side node (0, .5, .5) also lies on the slanted
face, whose share cancels its own. So the expected summed reaction is
(2/3) * 0.5 * P. The x = 0 face is held by a prescribed zero displacement,
because FEBio logs no reactions on "zero displacement" dofs. Lateral nodes
must contract to lam^-1/2. The material is
the writer's neo-Hookean mapping (uncoupled Mooney-Rivlin, c1 = mu/2, c2 = 0)
with k = 1000 mu, so near-incompressibility costs well under 1% here.

Skipped unless febio_exe is configured (simulation/local.toml or FEBIO_EXE).
"""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from simulation.core.materials import Fit
from simulation.stage1_finger.backends.febio import febio_exe, febio_material, run_febio
from simulation.stage1_finger.backends.febio.results import log_status, read_node_data

pytestmark = pytest.mark.skipif(febio_exe() is None, reason="febio_exe not configured")

MU_KPA = 40.0
LAM = 1.5
CORNERS = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]


def _nodes():
    return np.vstack([CORNERS] + [(CORNERS[i] + CORNERS[j]) / 2 for i, j in EDGES])


def _num(v):
    return f"{float(v):.12g}"


def _sub(parent, tag, text=None, **attrib):
    el = ET.SubElement(parent, tag, attrib)
    if text is not None:
        el.text = str(text)
    return el


def _feb(path):
    x = _nodes()
    mtype, params, _ = febio_material(Fit("neo_hookean", {"mu_kPa": MU_KPA}, "", "low"), 1000.0)
    root = ET.Element("febio_spec", version="4.0")
    _sub(root, "Module", type="solid")
    c = _sub(root, "Control")
    _sub(c, "analysis", "STATIC")
    _sub(c, "time_steps", 10)
    _sub(c, "step_size", 0.1)
    ts = _sub(c, "time_stepper", type="default")
    _sub(ts, "dtmax", 0.1)
    s = _sub(c, "solver", type="solid")
    _sub(s, "max_refs", 25)
    mat = _sub(_sub(root, "Material"), "material", id="1", name="m", type=mtype)
    for k, v in params.items():
        _sub(mat, k, _num(v))
    mesh = _sub(root, "Mesh")
    nodes = _sub(mesh, "Nodes", name="all")
    for i, p in enumerate(x, start=1):
        _sub(nodes, "node", ",".join(map(_num, p)), id=str(i))
    _sub(_sub(mesh, "Elements", type="tet10", name="e"), "elem", ",".join(map(str, range(1, 11))),
         id="1")
    ids = lambda mask: ",".join(str(i + 1) for i in np.flatnonzero(mask))
    _sub(mesh, "NodeSet", ids(x[:, 0] == 0), name="x0")
    _sub(mesh, "NodeSet", ids(x[:, 1] == 0), name="y0")
    _sub(mesh, "NodeSet", ids(x[:, 2] == 0), name="z0")
    _sub(_sub(root, "MeshDomains"), "SolidDomain", name="e", mat="m")
    b = _sub(root, "Boundary")
    bc = _sub(b, "bc", name="x0", node_set="x0", type="prescribed displacement")
    _sub(bc, "dof", "x")
    _sub(bc, "value", "0", lc="1")
    _sub(bc, "relative", 0)
    for name, dof in (("y0", "y"), ("z0", "z")):
        bc = _sub(b, "bc", name=name, node_set=name, type="zero displacement")
        for d in "xyz":
            _sub(bc, f"{d}_dof", int(d == dof))
    # every pulled node gets ux = (lam - 1) x, written as one bc per distinct x
    for xv in sorted(set(x[x[:, 0] > 0, 0])):
        name = f"pull_{xv:g}".replace(".", "p")
        _sub(mesh, "NodeSet", ids(x[:, 0] == xv), name=name)
        bc = _sub(b, "bc", name=name, node_set=name, type="prescribed displacement")
        _sub(bc, "dof", "x")
        _sub(bc, "value", _num((LAM - 1) * xv), lc="1")
        _sub(bc, "relative", 0)
    lc = _sub(_sub(root, "LoadData"), "load_controller", id="1", name="ramp", type="loadcurve")
    _sub(lc, "interpolate", "LINEAR")
    pts = _sub(lc, "points")
    _sub(pts, "pt", "0,0")
    _sub(pts, "pt", "1,1")
    log = _sub(_sub(root, "Output"), "logfile")
    _sub(log, "node_data", name="x0", node_set="x0", file="x0.txt", delim=",", data="Rx;Ry;Rz")
    _sub(mesh, "NodeSet", "3", name="ycorner")
    _sub(log, "node_data", name="ycorner", node_set="ycorner", file="y1.txt", delim=",",
         data="uy")
    ET.indent(root)
    path.write_text('<?xml version="1.0" encoding="ISO-8859-1"?>\n'
                    + ET.tostring(root, encoding="unicode").replace(" />", "/>"), encoding="latin-1")


def test_uniaxial_neo_hookean(tmp_path):
    feb = tmp_path / "uniaxial.feb"
    _feb(feb)
    info = run_febio(feb, timeout_s=120)
    assert info["returncode"] == 0, (tmp_path / "febio_stdout.txt").read_text()
    assert log_status(tmp_path / "uniaxial.log") == "normal"
    last = read_node_data(tmp_path / "x0.txt")[-1]
    assert last["time"] == pytest.approx(1.0)
    # FEBio logs the force the body exerts on the support: the x = 0 face of a
    # bar pulled along +x reads +x. sweep.py's blocked-force sign relies on this.
    force_N = last["values"][:, 0].sum()
    assert force_N > 0
    P_MPa = MU_KPA * 1e-3 * (LAM - LAM ** -2)
    assert force_N == pytest.approx(2 / 3 * 0.5 * P_MPa, rel=0.01)
    uy = read_node_data(tmp_path / "y1.txt")[-1]["values"][0, 0]
    assert 1 + uy == pytest.approx(LAM ** -0.5, rel=1e-3)
