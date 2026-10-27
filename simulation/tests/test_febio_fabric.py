"""Fabric layer material (sheet card in a silicone matrix), solved by FEBio on one tet10.

The material is exactly what the writer emits for a "cloth" part: the
fabric_woven_cotton card cast into Ecoflex 00-30 (fit 0) as a coupled solid
mixture of a neo-Hookean matrix and tension-only fiber-pow-linear warp (x)
and weft (y) families, for a 1 mm modelled layer (backends/febio/materials.py).

Same element and boundary scheme as test_febio_uniaxial.py: the corner tet
(0,0,0) (1,0,0) (0,1,0) (0,0,1) mm, symmetry planes on the three coordinate
faces, every node off the pulled axis' zero plane given u_a = (lam - 1) x_a,
and the other two directions free. The state is homogeneous uniaxial stress,
which tet10 represents exactly, and the summed reaction on the held face is
(2/3) * 0.5 mm^2 * P, with P the nominal stress along the pull (see that file).

Analytic expectation, principal axes aligned with the fibres, F = diag(l1, l2, l3):
  sigma_ii = mu/J (l_i^2 - 1) + lame ln(J)/J + tau(l_i)/J   (FEBio coupled neo-Hookean
             plus a fibre term on the warp and weft axes when l_i > 1)
  tau(l)   = l * P_f(l), P_f the fibre's nominal stress (materials.fiber_pow_linear_nominal)
and the two lateral stretches solve sigma = 0 on the free faces.

Checks: tension along x (warp) and y (weft) in the linear region and in the
toe, against that expectation; compression along x and y, where the pulled
fibre is slack and only the matrix (plus the other, now stretched, fibre
family) resists; and the physics: a layer per-width tension stiffness that
matches the card, and a tension/compression stiffness ratio in the thousands.

Skipped unless febio_exe is configured.
"""

import xml.etree.ElementTree as ET

import numpy as np
import pytest
from scipy.optimize import fsolve

from simulation.core.materials import load_material
from simulation.stage1_finger.backends.febio import febio_exe, run_febio
from simulation.stage1_finger.backends.febio.materials import (
    febio_sheet_mixture, fiber_pow_linear_c, fiber_pow_linear_nominal)
from simulation.stage1_finger.backends.febio.results import log_status, read_node_data

T_MM = 1.0
CORNERS = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
X = np.vstack([CORNERS] + [(CORNERS[i] + CORNERS[j]) / 2 for i, j in EDGES])


def _mixture():
    sheet = load_material("fabric_woven_cotton").fits[0]
    matrix = load_material("ecoflex_00-30").fits[0]
    return febio_sheet_mixture(sheet, matrix, T_MM, "x", "y", 1.0, 100.0), sheet, matrix


def _analytic(axis, lam):
    """(P along the pull in MPa, stretches (l1, l2, l3)) for uniaxial stress."""
    (_, parts, _), _, _ = _mixture()
    E, v = parts[0][1]["E"], parts[0][1]["v"]
    mu, lame = E / (2 * (1 + v)), v * E / ((1 + v) * (1 - 2 * v))
    fib = {0: parts[1][1], 1: parts[2][1]}      # warp on x, weft on y

    def sigma(l):
        J = np.prod(l)
        s = mu / J * (l ** 2 - 1) + lame * np.log(J) / J
        for i, f in fib.items():
            s[i] += l[i] * fiber_pow_linear_nominal(l[i], f["E"], f["lam0"], f["beta"]) / J
        return s

    lat = [i for i in range(3) if i != axis]

    def resid(q):
        l = np.ones(3)
        l[axis], l[lat[0]], l[lat[1]] = lam, q[0], q[1]
        return sigma(l)[lat] / mu

    q = fsolve(resid, [lam ** -0.5] * 2, xtol=1e-13)
    l = np.ones(3)
    l[axis], l[lat[0]], l[lat[1]] = lam, q[0], q[1]
    assert np.abs(resid(q)).max() < 1e-8
    return np.prod(l) * sigma(l)[axis] / lam, l


def _num(v):
    return f"{float(v):.12g}"


def _sub(parent, tag, text=None, **attrib):
    el = ET.SubElement(parent, tag, attrib)
    if text is not None:
        el.text = str(text)
    return el


def _feb(path, axis, lam):
    (mtype, parts, _), _, matrix = _mixture()
    root = ET.Element("febio_spec", version="4.0")
    _sub(root, "Module", type="solid")
    c = _sub(root, "Control")
    _sub(c, "analysis", "STATIC")
    _sub(c, "time_steps", 20)
    _sub(c, "step_size", 0.05)
    _sub(_sub(c, "time_stepper", type="default"), "dtmax", 0.05)
    _sub(_sub(c, "solver", type="solid"), "max_refs", 25)
    m = _sub(_sub(root, "Material"), "material", id="1", name="m", type=mtype)
    _sub(m, "density", "1.07e-9")
    for ctype, params in parts:
        s = _sub(m, "solid", type=ctype)
        for k, v in params.items():
            if k == "fiber":
                _sub(s, "fiber", v, type="vector")
            else:
                _sub(s, k, _num(v))
    mesh = _sub(root, "Mesh")
    nodes = _sub(mesh, "Nodes", name="all")
    for i, p in enumerate(X, start=1):
        _sub(nodes, "node", ",".join(map(_num, p)), id=str(i))
    _sub(_sub(mesh, "Elements", type="tet10", name="e"), "elem",
         ",".join(map(str, range(1, 11))), id="1")
    ids = lambda mask: ",".join(str(i + 1) for i in np.flatnonzero(mask))
    d = "xyz"
    _sub(_sub(root, "MeshDomains"), "SolidDomain", name="e", mat="m")
    b = _sub(root, "Boundary")
    for k in range(3):
        _sub(mesh, "NodeSet", ids(X[:, k] == 0), name=f"{d[k]}0")
        if k == axis:     # the held face: prescribed zero, so FEBio logs its reactions
            bc = _sub(b, "bc", name=f"{d[k]}0", node_set=f"{d[k]}0", type="prescribed displacement")
            _sub(bc, "dof", d[k])
            _sub(bc, "value", "0", lc="1")
            _sub(bc, "relative", 0)
        else:
            bc = _sub(b, "bc", name=f"{d[k]}0", node_set=f"{d[k]}0", type="zero displacement")
            for dd in d:
                _sub(bc, f"{dd}_dof", int(dd == d[k]))
    for xv in sorted(set(X[X[:, axis] > 0, axis])):
        name = f"pull_{xv:g}".replace(".", "p")
        _sub(mesh, "NodeSet", ids(X[:, axis] == xv), name=name)
        bc = _sub(b, "bc", name=name, node_set=name, type="prescribed displacement")
        _sub(bc, "dof", d[axis])
        _sub(bc, "value", _num((lam - 1) * xv), lc="1")
        _sub(bc, "relative", 0)
    lc = _sub(_sub(root, "LoadData"), "load_controller", id="1", name="ramp", type="loadcurve")
    _sub(lc, "interpolate", "LINEAR")
    pts = _sub(lc, "points")
    _sub(pts, "pt", "0,0")
    _sub(pts, "pt", "1,1")
    log = _sub(_sub(root, "Output"), "logfile")
    _sub(log, "node_data", name="held", node_set=f"{d[axis]}0", file="held.txt", delim=",",
         data="Rx;Ry;Rz")
    for k, node in ((0, 2), (1, 3), (2, 4)):     # corner nodes (1,0,0), (0,1,0), (0,0,1)
        _sub(mesh, "NodeSet", str(node), name=f"c{d[k]}")
        _sub(log, "node_data", name=f"c{d[k]}", node_set=f"c{d[k]}", file=f"c{d[k]}.txt",
             delim=",", data=f"u{d[k]}")
    ET.indent(root)
    path.write_text('<?xml version="1.0" encoding="ISO-8859-1"?>\n'
                    + ET.tostring(root, encoding="unicode").replace(" />", "/>"), encoding="latin-1")


def _solve(tmp_path, axis, lam):
    feb = tmp_path / "fabric.feb"
    _feb(feb, axis, lam)
    info = run_febio(feb, timeout_s=120)
    assert info["returncode"] == 0, (tmp_path / "febio_stdout.txt").read_text()
    assert log_status(tmp_path / "fabric.log") == "normal"
    last = read_node_data(tmp_path / "held.txt")[-1]
    assert last["time"] == pytest.approx(1.0)
    force = last["values"][:, axis].sum()
    stretch = [1 + read_node_data(tmp_path / f"c{d}.txt")[-1]["values"][0, 0] for d in "xyz"]
    return force, np.array(stretch)


def test_mixture_parameters():
    (mtype, parts, info), sheet, matrix = _mixture()
    assert mtype == "solid mixture"
    assert [p[0] for p in parts] == ["neo-Hookean", "fiber-pow-linear", "fiber-pow-linear"]
    assert parts[1][1]["fiber"] == "1,0,0" and parts[2][1]["fiber"] == "0,1,0"
    # matrix keeps the silicone's mu0 and K
    E, v = parts[0][1]["E"], parts[0][1]["v"]
    assert E / (2 * (1 + v)) == pytest.approx(matrix.initial_shear_modulus_kPa() * 1e-3)
    assert E / (3 * (1 - 2 * v)) == pytest.approx(100 * matrix.initial_shear_modulus_kPa() * 1e-3)
    # fibre E chosen so the layer's linear-region dT/dlam equals the card's S
    for k, S in zip((1, 2), sheet.warp_weft_N_per_mm()):
        f = parts[k][1]
        assert T_MM * f["E"] * fiber_pow_linear_c(f["lam0"], f["beta"]) == pytest.approx(S)
        lam = np.array([1.02, 1.03])
        P = fiber_pow_linear_nominal(lam, f["E"], f["lam0"], f["beta"])
        assert T_MM * np.diff(P)[0] / 0.01 == pytest.approx(S, rel=1e-9)
    # tension only, continuous at lam0
    f = parts[1][1]
    assert fiber_pow_linear_nominal(0.95, f["E"], f["lam0"], f["beta"]) == 0.0
    lo, hi = fiber_pow_linear_nominal([f["lam0"] - 1e-9, f["lam0"] + 1e-9], f["E"], f["lam0"],
                                      f["beta"])
    assert lo == pytest.approx(hi, rel=1e-6)
    # scale and thickness: per-width stiffness scales, fibre modulus ~ 1/t
    sheet_fit = load_material("fabric_woven_cotton").fits[0]
    _, p2, i2 = febio_sheet_mixture(sheet_fit, matrix, 2 * T_MM, "x", "y", 4.0, 100.0)
    assert p2[1][1]["E"] == pytest.approx(2 * parts[1][1]["E"])
    assert i2["warp_N_per_mm"] == pytest.approx(4 * sheet.warp_weft_N_per_mm()[0])


febio = pytest.mark.skipif(febio_exe() is None, reason="febio_exe not configured")

CASES = [
    # (pull axis, stretch, regime)
    (0, 1.02, "warp linear"),
    (1, 1.02, "weft linear"),
    (0, 1.003, "warp toe"),
    (0, 0.90, "x compression"),
    (1, 0.90, "y compression"),
]


@febio
@pytest.mark.parametrize("axis,lam,regime", CASES, ids=[c[2] for c in CASES])
def test_single_element_against_analytic(tmp_path, axis, lam, regime):
    force, stretch = _solve(tmp_path, axis, lam)
    P, l = _analytic(axis, lam)
    assert force == pytest.approx(2 / 3 * 0.5 * P, rel=0.01)
    assert stretch == pytest.approx(l, abs=2e-4)
    if lam < 1:
        assert force < 0      # compression pushes on the held face


@febio
def test_stiff_in_tension_soft_in_compression(tmp_path):
    """Per-width tension stiffness matches the card; compression is matrix-only."""
    sheet = load_material("fabric_woven_cotton").fits[0]
    S_warp, S_weft = sheet.warp_weft_N_per_mm()
    out = {}
    for name, axis, lam in (("x+", 0, 1.02), ("x++", 0, 1.03), ("y+", 1, 1.02), ("y++", 1, 1.03),
                            ("x-", 0, 0.98)):
        d = tmp_path / name.replace("+", "p").replace("-", "m")
        d.mkdir()
        f, _ = _solve(d, axis, lam)
        out[name] = f / (2 / 3 * 0.5)            # nominal stress P in MPa
    # slope between 2% and 3% stretch, times the 1 mm layer, is the per-width stiffness
    # (the matrix adds ~0.1 N/mm)
    assert (out["x++"] - out["x+"]) / 0.01 * T_MM == pytest.approx(S_warp, rel=0.02)
    assert (out["y++"] - out["y+"]) / 0.01 * T_MM == pytest.approx(S_weft, rel=0.02)
    k_tension = out["x+"] / 0.02
    k_compression = out["x-"] / -0.02
    assert 0 < k_compression < 0.2                # MPa: plane-strain-like matrix, ~4 mu0
    assert k_tension / k_compression > 50
