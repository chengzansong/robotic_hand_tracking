"""Canonical material fits (simulation/materials/README.md) -> FEBio 4 materials.

All targets are FEBio's uncoupled (deviatoric + volumetric) hyperelastic
materials, which suit near-incompressible rubber with tet10 elements. Units
in: kPa. Units out: MPa (solver set mm, N, MPa).

Mappings, checked against the FEBio 4 user manual (help.febio.org/docs/FEBioUser-4-7/,
sections 4.1.2.9 Mooney-Rivlin, 4.1.2.11 Ogden, 4.1.2.20 Yeoh) and the v4.13
source (FEBioMech/FEMooneyRivlin.cpp, FEOgdenMaterial.cpp, FEYeoh.cpp):

  card model    FEBio type       FEBio W (deviatoric)                 conversion
  neo_hookean   "Mooney-Rivlin"  c1 (I1~-3) + c2 (I2~-3)              c1 = mu/2, c2 = 0
  mooney_rivlin "Mooney-Rivlin"  same                                 c1 = C10, c2 = C01
  yeoh          "Yeoh"           sum_i c_i (I1~-3)^i, i = 1..6        c_i = C_i0
  ogden         "Ogden"          sum c_i/m_i^2 (sum_k lam~_k^m_i - 3) c_i = 2 mu_i, m_i = alpha_i

Volumetric part (FEUncoupledMaterial, pressure_model "default"):
U = k/2 (ln J)^2. The bulk modulus k comes from the fit's bulk_modulus_kPa,
else from its poisson ratio, else from the study's bulk_to_shear ratio.
"""

KPA_TO_MPA = 1e-3
MAX_TERMS = 6   # FEBio Ogden and Yeoh take c1..c6
SUPPORTED = ("neo_hookean", "mooney_rivlin", "yeoh", "ogden")


class UnsupportedMaterial(ValueError):
    pass


def bulk_modulus_kPa(fit, bulk_to_shear):
    """Bulk modulus for the volumetric term, and where it came from."""
    mu0 = fit.initial_shear_modulus_kPa()
    if fit.bulk_modulus_kPa is not None:
        return float(fit.bulk_modulus_kPa), "card bulk_modulus_kPa"
    if fit.poisson is not None:
        nu = float(fit.poisson)
        if not 0 < nu < 0.5:
            raise ValueError(f"poisson {nu} must be in (0, 0.5)")
        return 2 * mu0 * (1 + nu) / (3 * (1 - 2 * nu)), f"card poisson {nu}"
    return bulk_to_shear * mu0, f"study bulk_to_shear {bulk_to_shear}"


def febio_material(fit, bulk_to_shear):
    """(FEBio type string, {param: value in MPa or dimensionless}, k source)."""
    p = fit.params
    if fit.model not in SUPPORTED:
        raise UnsupportedMaterial(
            f"no FEBio mapping for model {fit.model!r}; supported: {', '.join(SUPPORTED)}")
    k_kPa, k_source = bulk_modulus_kPa(fit, bulk_to_shear)
    if fit.model == "neo_hookean":
        mtype, params = "Mooney-Rivlin", {"c1": p["mu_kPa"] / 2, "c2": 0.0}
    elif fit.model == "mooney_rivlin":
        mtype, params = "Mooney-Rivlin", {"c1": p["c10_kPa"], "c2": p["c01_kPa"]}
    elif fit.model == "yeoh":
        mtype = "Yeoh"
        params = {"c1": p["c10_kPa"], "c2": p["c20_kPa"], "c3": p["c30_kPa"]}
    elif fit.model == "ogden":
        mu, alpha = list(p["mu_kPa"]), list(p["alpha"])
        if len(mu) > MAX_TERMS:
            raise UnsupportedMaterial(f"FEBio Ogden takes at most {MAX_TERMS} terms, got {len(mu)}")
        if any(a == 0 for a in alpha):
            raise UnsupportedMaterial("FEBio Ogden requires every alpha (m_i) to be nonzero")
        mtype, params = "Ogden", {}
        for i, (m, a) in enumerate(zip(mu, alpha), start=1):
            params[f"c{i}"] = 2 * m
            params[f"m{i}"] = a

    out = {}
    for name, value in params.items():
        out[name] = value if name.startswith("m") else value * KPA_TO_MPA
    out["k"] = k_kPa * KPA_TO_MPA
    return mtype, out, k_source


def febio_initial_shear_MPa(mtype, params):
    """Small-strain shear modulus of an FEBio uncoupled material."""
    if mtype == "Mooney-Rivlin":
        return 2 * (params["c1"] + params["c2"])
    if mtype == "Yeoh":
        return 2 * params["c1"]
    if mtype == "Ogden":
        # d2W gives mu0 = sum c_i / 2 for W = sum c_i/m_i^2 (...)
        return sum(v for k, v in params.items() if k.startswith("c")) / 2
    raise UnsupportedMaterial(mtype)
