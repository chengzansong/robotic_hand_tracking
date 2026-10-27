"""Write an FEBio 4.0 input file (.feb) from an FEMesh, material fits, and a load case.

Syntax follows the FEBio 4 user manual (help.febio.org/docs/FEBioUser-4-7/,
chapter 3) where it matches the v4.13 source and the 4.0-format regression
files in github.com/febiosoftware/TestSuite (Verify/ma06.feb is the skeleton
used here). Where the manual and source disagree, the source wins; those
places are commented.

Model:
- Module "solid", STATIC analysis, auto time stepping over t in [0, 1].
- One tet10 Elements block and SolidDomain per part. With an uncoupled
  material and tet10, FEBio 4.13 picks the standard displacement
  "elastic-solid" domain (FESolidDomainFactory; FEModelBuilder sets
  m_b3field_tet = false). Quadratic tets do not lock the way tet4 does, so
  no three-field formulation is needed. [solver] domain_type overrides it.
- A sheet material (fabric card) on a layer part is written as a coupled
  "solid mixture": neo-Hookean matrix plus tension-only fiber-pow-linear warp
  and weft families along global axes (materials.febio_sheet_mixture). Its
  domain is the same tet10 elastic-solid; parts are bonded by shared nodes.
- Pressure on the cavity surface, ramped linearly by load controller 1.
  Facets point out of the solid, so a positive pressure inflates
  (FEPressureLoad: traction = -p * N).
- Clamp node set held in x, y, z by prescribed zero displacements, so its
  reactions are logged (FEBio reports none for "zero displacement").
- blocked_force case: the tip node set is also held by a prescribed zero
  displacement along the bending axis, and its reaction is the blocked force.
  This holds the whole tip face, so the held nodes can also apply a couple
  once the face rotates (1.43 N mm at 5 kPa on SW-6, mesh A).
- blocked_line case: a line support instead, like a fingertip on a load cell.
  Only the tip-face nodes on the bend-side edge of the face (the underside
  corner line across the width, see tip_line_nodes) are held along the bending
  axis, as node set "tip_line"; the tip is free to rotate about that line and
  to slide along the finger. The tip log still covers the whole face. A free
  dof logs its internal-force residual, which is about zero where no load acts
  (checked on a small beam in tests/test_blocked_line.py). No pressure or
  contact touches the tip face, so the summed tip reaction is the support force.
- Half model (mesh symmetry "half_y"): uy = 0 on the symmetry node set.
- Contact between the facing walls of each slot between chambers: one
  SurfacePair and one "sliding-elastic" interface per slot, walls wound out
  of the solid, two-pass, auto penalty, frictionless unless fric_coeff is set.
  Parameters are those of FEBioMech/FESlidingElasticInterface.cpp (v4.13)
  and match TestSuite Verify/co*.feb. Contact pressure and gap go to the plot
  file. Tip curling onto the base is not modelled.
- Must points at each requested pressure level (load controller 2, used as
  the dtmax curve), so output lands exactly on them. The curve is STEP
  interpolated and PointCurve::value returns the next point's value
  (FECore/PointCurve.cpp), so point (t_k, t_k - t_k-1) caps each interval
  at its own length.
- Log files (written at must points only): tip node positions, displacements,
  and reactions; clamp reactions; and displacements of every node, all
  comma-delimited text.
"""

import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from .materials import febio_material, febio_sheet_mixture

KPA_TO_MPA = 1e-3
KG_M3_TO_TONNE_MM3 = 1e-12

# Defaults are the settings tuned on the smoke mesh (README, "Solver tuning"):
# full Newton with a non-symmetric tangent, 2 initial steps per output level,
# aggressive step growth. Full Newton cut wall time 20-35 % against BFGS; the
# non-symmetric tangent is what lets sliding-elastic contact converge (the
# symmetric run stalled at 2.35 kPa on max reformations).
SOLVER_DEFAULTS = {
    "substeps": 2,            # initial time steps per output level
    "max_retries": 10,
    "opt_iter": 15,
    "aggressiveness": 1,      # time stepper: 1 grows the step faster after easy steps
    "max_refs": 25,
    "max_ups": 0,             # BFGS updates; 0 = full Newton (FECore/BFGSSolver.cpp)
    "symmetric_stiffness": "non-symmetric",   # or "symmetric" (faster without contact)
    "pressure_symmetric": False,          # symmetric_stiffness flag of the pressure load
    "threads": 24,            # OMP_NUM_THREADS for the solver; 0 = leave unset
    "dtol": 0.001,
    "etol": 0.01,
    "rtol": 0.0,
    "lstol": 0.9,
    "domain_type": "",        # "" lets FEBio choose (elastic-solid for uncoupled tet10)
    "elem_type": "",          # "" = FEBio default TET10G8; e.g. "TET10G4"
}

CONTACT_DEFAULTS = {
    "enabled": True,
    "laugon": "PENALTY",      # PENALTY or AUGLAG
    "penalty": 1.0,           # scale factor on the automatic penalty
    "auto_penalty": True,
    "two_pass": True,
    "tolerance": 0.2,         # augmentation tolerance (AUGLAG only)
    "search_tol": 0.01,
    "search_radius": 1.0,     # mm; also the deepest penetration still treated as contact
    "symmetric_stiffness": False,
    "fric_coeff": 0.0,        # silicone is sticky; frictionless is the first pass
    "penetration_tol_mm": 0.05,   # post-solve check only, not written to FEBio
}

LOG_FILES = {
    "tip": "tip.txt",
    "clamp": "clamp.txt",
    "all_nodes": "nodes_u.txt",
}


def _fmt(v):
    return f"{v:.10g}"


def _sub(parent, tag, text=None, **attrib):
    el = ET.SubElement(parent, tag, {k: str(v) for k, v in attrib.items()})
    if text is not None:
        el.text = text if isinstance(text, str) else _fmt(text)
    return el


def _id_list(ids):
    """Comma-separated 1-based ids, wrapped so lines stay short. FEBio 4.0 needs
    a text list in <NodeSet> (child <n> tags throw in FEBioMeshSection4.cpp)."""
    ids = [str(int(i) + 1) for i in ids]
    return "\n" + ",\n".join(", ".join(ids[i:i + 16]) for i in range(0, len(ids), 16)) + "\n"


def _control(root, fractions, opts, end_time=1.0):
    dt0 = fractions[0] / opts["substeps"]
    n = max(1, round(end_time / dt0))
    c = _sub(root, "Control")
    _sub(c, "analysis", "STATIC")
    _sub(c, "time_steps", str(n))
    _sub(c, "step_size", end_time / n)
    _sub(c, "plot_zero_state", "1")
    _sub(c, "plot_level", "PLOT_MUST_POINTS")
    # Logfile records only at must points, so the all-node log stays one record
    # per requested pressure level (enum in FECore/FEAnalysis.cpp, used by
    # FEBioModel::Write for the data store, v4.13).
    _sub(c, "output_level", "OUTPUT_MUST_POINTS")
    ts = _sub(c, "time_stepper", type="default")
    _sub(ts, "max_retries", str(opts["max_retries"]))
    _sub(ts, "opt_iter", str(opts["opt_iter"]))
    _sub(ts, "dtmin", dt0 / 1000)
    # In 4.x a load curve on dtmax makes its points must points; the curve
    # value is the max step. A nonzero text value is ignored with a warning
    # (FECore/FETimeStepController.cpp), so write 0 as the TestSuite does.
    _sub(ts, "dtmax", "0", lc="2")
    _sub(ts, "aggressiveness", str(opts["aggressiveness"]))
    _sub(ts, "cutback", "0.5")
    s = _sub(c, "solver", type="solid")
    _sub(s, "symmetric_stiffness", opts["symmetric_stiffness"])
    _sub(s, "lstol", opts["lstol"])
    _sub(s, "max_refs", str(opts["max_refs"]))
    _sub(s, "reform_each_time_step", "1")
    _sub(s, "diverge_reform", "1")
    _sub(s, "dtol", opts["dtol"])
    _sub(s, "etol", opts["etol"])
    _sub(s, "rtol", opts["rtol"])
    # max_ups belongs to the quasi-Newton method in 4.0 (FECore/BFGSSolver.cpp),
    # not directly under <solver> as the manual's example shows.
    qn = _sub(s, "qn_method", type="BFGS")
    _sub(qn, "max_ups", str(opts["max_ups"]))


def _materials(root, part_materials, bulk_to_shear):
    mat_el = _sub(root, "Material")
    used = {}
    for i, (part, pm) in enumerate(part_materials.items(), start=1):
        name = f"{part}_material"
        if getattr(pm, "is_sheet", False):
            used[part] = _sheet_material(mat_el, i, name, pm, bulk_to_shear)
            continue
        mtype, params, k_source = febio_material(pm.fit, bulk_to_shear)
        m = _sub(mat_el, "material", id=str(i), name=name, type=mtype)
        _sub(m, "density", pm.material.density_kg_m3 * KG_M3_TO_TONNE_MM3)
        for key in sorted(params, key=lambda k: (k[0] != "c", k)):
            if key != "k":
                _sub(m, key, params[key])
        _sub(m, "k", params["k"])
        _sub(m, "pressure_model", "default")
        used[part] = {"name": name, "type": mtype, "params_MPa": params, "k_source": k_source}
    return used


def _sheet_material(mat_el, i, name, pm, bulk_to_shear):
    """A sheet (fabric) cast into its matrix: a coupled "solid mixture" of a
    neo-Hookean matrix and tension-only fiber-pow-linear warp and weft families
    (see materials.febio_sheet_mixture). The fibre vectors are global axes: the
    layer parts are flat, so no local material axes are needed."""
    layer = pm.layer or {}
    mtype, parts, info = febio_sheet_mixture(
        pm.fit, pm.matrix.fit, layer.get("thickness_mm"), layer.get("warp_axis", "x"),
        layer.get("weft_axis", "y"), pm.stiffness_scale, bulk_to_shear)
    m = _sub(mat_el, "material", id=str(i), name=name, type=mtype)
    _sub(m, "density", pm.matrix.material.density_kg_m3 * KG_M3_TO_TONNE_MM3)
    for ctype, params in parts:
        s = _sub(m, "solid", type=ctype)
        for key, value in params.items():
            if key == "fiber":
                _sub(s, "fiber", value, type="vector")
            else:
                _sub(s, key, value)
    return {"name": name, "type": mtype, "constituents": [{"type": t, **p} for t, p in parts],
            "params_MPa": info, "k_source": info["matrix_k_source"]}


def _tri6_surface(parent, name, facets):
    surf = _sub(parent, "Surface", name=name)
    for i, f in enumerate(facets, start=1):
        _sub(surf, "tri6", ",".join(str(int(n) + 1) for n in f), id=str(i))


def tip_line_nodes(mesh, bend_toward="-z", tol_mm=None):
    """Node ids of the blocked_line support: the tip-face nodes on the face's
    edge toward the bend side, i.e. at the minimum of z on the tip face for
    bend_toward = "-z" (the underside corner line, across the finger width).
    tol_mm defaults to the mesh's fem_tags plane_tol_mm."""
    if tol_mm is None:
        tol_mm = (mesh.info or {}).get("fem_tags", {}).get("plane_tol_mm", 0.01)
    tip = np.asarray(mesh.node_sets["tip"])
    h = mesh.nodes[tip, "xyz".index(bend_toward.lstrip("+-"))]
    edge = h.max() if bend_toward.startswith("+") else h.min()
    return tip[np.abs(h - edge) <= tol_mm]


def _mesh(root, mesh, contact, extra_sets=None):
    m = _sub(root, "Mesh")
    nodes = _sub(m, "Nodes", name="all")
    for i, x in enumerate(mesh.nodes, start=1):
        _sub(nodes, "node", ",".join(_fmt(v) for v in x), id=str(i))
    eid = 1
    for part, tets in mesh.parts.items():
        el = _sub(m, "Elements", type="tet10", name=part)
        for t in tets:
            _sub(el, "elem", ",".join(str(int(n) + 1) for n in t), id=str(eid))
            eid += 1
    for name in ("clamp", "tip", "symmetry"):
        if name in mesh.node_sets:
            _sub(m, "NodeSet", _id_list(mesh.node_sets[name]), name=name)
    for name, ids in (extra_sets or {}).items():
        _sub(m, "NodeSet", _id_list(ids), name=name)
    _sub(m, "NodeSet", _id_list(np.arange(len(mesh.nodes))), name="all_nodes")
    _tri6_surface(m, "cavity", mesh.surfaces["cavity"])
    for slot in _contact_slots(mesh, contact):
        _tri6_surface(m, f"{slot}_lo", mesh.surfaces[f"{slot}_lo"])
        _tri6_surface(m, f"{slot}_hi", mesh.surfaces[f"{slot}_hi"])
        pair = _sub(m, "SurfacePair", name=slot)
        _sub(pair, "primary", f"{slot}_lo")
        _sub(pair, "secondary", f"{slot}_hi")


def _contact_slots(mesh, contact):
    if not contact["enabled"]:
        return []
    return [s for s in mesh.info.get("slots", [])
            if len(mesh.surfaces.get(f"{s}_lo", [])) and len(mesh.surfaces.get(f"{s}_hi", []))]


def _contact(root, mesh, contact):
    slots = _contact_slots(mesh, contact)
    if not slots:
        return
    sec = _sub(root, "Contact")
    for slot in slots:
        c = _sub(sec, "contact", name=slot, surface_pair=slot, type="sliding-elastic")
        _sub(c, "laugon", contact["laugon"])
        _sub(c, "tolerance", contact["tolerance"])
        _sub(c, "penalty", contact["penalty"])
        _sub(c, "auto_penalty", str(int(contact["auto_penalty"])))
        _sub(c, "two_pass", str(int(contact["two_pass"])))
        _sub(c, "search_tol", contact["search_tol"])
        _sub(c, "search_radius", contact["search_radius"])
        _sub(c, "symmetric_stiffness", str(int(contact["symmetric_stiffness"])))
        _sub(c, "tension", "0")
        _sub(c, "fric_coeff", contact["fric_coeff"])


def _domains(root, mesh, opts):
    d = _sub(root, "MeshDomains")
    for part in mesh.parts:
        attrib = {"name": part, "mat": f"{part}_material"}
        if opts["domain_type"]:
            attrib["type"] = opts["domain_type"]
        if opts["elem_type"]:
            # 4.0 strings are TET10G4, TET10G8, TET10GL11 (FEModelBuilder ElementSpec),
            # not the GAUSS4/GAUSS8 names in the manual.
            attrib["elem_type"] = opts["elem_type"]
        _sub(d, "SolidDomain", **attrib)


def _hold(b, name, node_set, dof):
    """Hold one dof at zero with a prescribed displacement of 0.

    FEBio 4.13 logs Rx/Ry/Rz as zero on nodes fixed by "zero displacement";
    only prescribed dofs report reactions (checked 2026-10-20 on a single
    tet10). Holding with a prescribed 0 gives the same constraint and a
    usable reaction force.
    """
    bc = _sub(b, "bc", name=name, node_set=node_set, type="prescribed displacement")
    _sub(bc, "dof", dof)
    _sub(bc, "value", "0", lc="1")
    _sub(bc, "relative", "0")


def _boundary(root, mesh, load, bend_dof):
    b = _sub(root, "Boundary")
    for dof in "xyz":
        _hold(b, f"clamp_{dof}", "clamp", dof)
    if "symmetry" in mesh.node_sets:
        bc = _sub(b, "bc", name="symmetry", node_set="symmetry", type="zero displacement")
        for dof in "xyz":
            _sub(bc, f"{dof}_dof", "1" if dof == "y" else "0")
    if load.case == "blocked_force":
        _hold(b, "tip_block", "tip", bend_dof)
    elif load.case == "blocked_line":
        _hold(b, "tip_line_block", "tip_line", bend_dof)


def _loads(root, load, opts):
    loads = _sub(root, "Loads")
    sl = _sub(loads, "surface_load", name="cavity_pressure", surface="cavity", type="pressure")
    _sub(sl, "pressure", load.max_kPa * KPA_TO_MPA, lc="1")
    _sub(sl, "linear", "0")
    _sub(sl, "symmetric_stiffness", str(int(bool(opts["pressure_symmetric"]))))


def _load_data(root, fractions):
    ld = _sub(root, "LoadData")
    lc = _sub(ld, "load_controller", id="1", name="pressure_ramp", type="loadcurve")
    _sub(lc, "interpolate", "LINEAR")
    _sub(lc, "extend", "CONSTANT")
    pts = _sub(lc, "points")
    _sub(pts, "pt", "0,0")
    _sub(pts, "pt", "1,1")
    mp = _sub(ld, "load_controller", id="2", name="must_points", type="loadcurve")
    _sub(mp, "interpolate", "STEP")
    _sub(mp, "extend", "CONSTANT")
    pts = _sub(mp, "points")
    prev = 0.0
    for t in fractions:
        _sub(pts, "pt", f"{_fmt(t)},{_fmt(t - prev)}")
        prev = t


def _output(root, with_contact):
    out = _sub(root, "Output")
    plot = _sub(out, "plotfile", type="febio")
    plot_vars = ["displacement", "stress", "relative volume", "reaction forces"]
    if with_contact:
        plot_vars += ["contact pressure", "contact gap"]
    for var in plot_vars:
        _sub(plot, "var", type=var)
    log = _sub(out, "logfile")
    _sub(log, "node_data", name="tip", node_set="tip", file=LOG_FILES["tip"], delim=",",
         data="x;y;z;ux;uy;uz;Rx;Ry;Rz")
    _sub(log, "node_data", name="clamp", node_set="clamp", file=LOG_FILES["clamp"],
         delim=",", data="Rx;Ry;Rz")
    _sub(log, "node_data", name="all_nodes", node_set="all_nodes",
         file=LOG_FILES["all_nodes"], delim=",", data="ux;uy;uz")


def _options(defaults, given, section):
    opts = dict(defaults)
    unknown = set(given or {}) - set(opts)
    if unknown:
        raise ValueError(f"unknown [{section}] keys for febio: {sorted(unknown)}")
    opts.update(given or {})
    return opts


def write_feb(path, mesh, part_materials, load, bulk_to_shear, solver=None, bend_dof="z",
              contact=None, end_time=1.0):
    """Write the .feb. Returns a dict describing what was written.

    end_time < 1 stops the analysis early at that load fraction, which is how
    the time estimate probes the first pressure level(s).
    """
    opts = _options(SOLVER_DEFAULTS, solver, "solver")
    contact = _options(CONTACT_DEFAULTS, contact, "contact")
    missing = set(mesh.parts) - set(part_materials)
    if missing:
        raise ValueError(f"no material for parts {sorted(missing)}")

    extra_sets = {}
    if load.case == "blocked_line":
        bend_toward = (mesh.info or {}).get("fem_tags", {}).get("bend_toward", f"-{bend_dof}")
        if bend_toward.lstrip("+-") != bend_dof:
            raise ValueError(f"bend_toward {bend_toward!r} is not along bend_dof {bend_dof!r}")
        extra_sets["tip_line"] = tip_line_nodes(mesh, bend_toward)

    root = ET.Element("febio_spec", version="4.0")
    root.append(ET.Comment(" Generated by simulation/stage1_finger/backends/febio/writer.py. "
                           "Units: mm, N, MPa, s, tonne/mm^3. "))
    _sub(root, "Module", type="solid")
    fractions = load.fractions()
    _control(root, fractions, opts, end_time)
    used = _materials(root, {p: part_materials[p] for p in mesh.parts}, bulk_to_shear)
    _mesh(root, mesh, contact, extra_sets)
    _domains(root, mesh, opts)
    _boundary(root, mesh, load, bend_dof)
    _loads(root, load, opts)
    _contact(root, mesh, contact)
    _load_data(root, fractions)
    slots = _contact_slots(mesh, contact)
    _output(root, bool(slots))

    ET.indent(root, space="  ")
    path = Path(path)
    # Same declaration and empty-tag style as files written by FEBio Studio.
    text = ET.tostring(root, encoding="unicode").replace(" />", "/>")
    path.write_text('<?xml version="1.0" encoding="ISO-8859-1"?>\n' + text + "\n",
                    encoding="latin-1")
    out = {"path": str(path), "materials": used, "solver": opts, "contact": contact,
           "contact_slots": slots, "log_files": LOG_FILES, "levels_kPa": load.levels(),
           "end_time": end_time}
    if "tip_line" in extra_sets:
        out["tip_line_nodes"] = len(extra_sets["tip_line"])
    return out
