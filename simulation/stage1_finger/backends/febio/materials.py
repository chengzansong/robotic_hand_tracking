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


# --- Sheet reinforcement (cast-in fabric) ------------------------------------
#
# A sheet card (kind = "sheet", model "fabric_pow_linear") on a layer part of
# modelled thickness t becomes a coupled FEBio "solid mixture" of
#   - the matrix silicone as a coupled "neo-Hookean" (E, v) with the matrix
#     fit's initial shear modulus mu0 and bulk modulus K (same K rule as the
#     uncoupled cards: card value, else Poisson, else bulk_to_shear * mu0);
#   - one "fiber-pow-linear" family along the warp axis and one along the weft
#     axis, tension only (FEBio skips a fibre with I_n < 1).
# Coupled, not uncoupled: an uncoupled fibre sees only the deviatoric stretch
# J^-1/3 lam, so with fibres much stiffer than the bulk modulus the layer could
# stretch by dilating, capped near 9 K whatever the fibre modulus. The coupled
# fibre uses the true stretch.
#
# FEBio 4.13 fiber-pow-linear (FEBioMech/FEFiberPowLinear.cpp), fibre invariant
# I = lam^2 along the reference direction n0, Kirchhoff fibre stress tau = J sigma_n:
#   toe    (1 <= I < I0): tau = 2 I ksi (I - 1)^(beta - 1)
#   linear (I >= I0):     tau = 2 b I - E sqrt(I)
#   ksi = E / (4 (beta - 1)) I0^(-3/2) (I0 - 1)^(2 - beta),  I0 = lam0^2
#   b   = ksi (I0 - 1)^(beta - 1) + E / (2 sqrt(I0))
# The nominal stress along the fibre, P = tau / lam, is then linear in lam past
# lam0 with slope dP/dlam = 2 b = E c, c = 1/lam0 + (lam0^2 - 1) / (2 (beta - 1) lam0^3).
# A layer of thickness t carries T = t P per unit width, so the card's
# linear-region stiffness per width S = dT/dlam gives E = S / (t c).

SHEET_MODELS = ("fabric_pow_linear",)
AXIS_VEC = {"x": (1.0, 0.0, 0.0), "y": (0.0, 1.0, 0.0), "z": (0.0, 0.0, 1.0)}


def fiber_pow_linear_c(lam0, beta):
    """dP/dlam in the linear region divided by FEBio's E."""
    return 1.0 / lam0 + (lam0 ** 2 - 1) / (2 * (beta - 1) * lam0 ** 3)


def fiber_pow_linear_nominal(lam, E, lam0, beta):
    """Nominal (1st Piola-Kirchhoff) fibre stress P(lam) of FEBio's
    fiber-pow-linear, in the units of E; zero for lam <= 1 (tension only)."""
    import numpy as np
    lam = np.asarray(lam, float)
    I, I0 = lam ** 2, lam0 ** 2
    ksi = E / 4 / (beta - 1) * I0 ** -1.5 * (I0 - 1) ** (2 - beta)
    b = ksi * (I0 - 1) ** (beta - 1) + E / 2 / np.sqrt(I0)
    tau = np.where(I < I0, 2 * I * ksi * np.clip(I - 1, 0, None) ** (beta - 1),
                   2 * b * I - E * np.sqrt(I))
    return np.where(I > 1, tau / lam, 0.0)


def coupled_neo_hookean(fit, bulk_to_shear):
    """(E, v) in MPa for FEBio's coupled "neo-Hookean" with the fit's mu0 and K."""
    mu = fit.initial_shear_modulus_kPa() * KPA_TO_MPA
    k_kPa, k_source = bulk_modulus_kPa(fit, bulk_to_shear)
    K = k_kPa * KPA_TO_MPA
    E = 9 * K * mu / (3 * K + mu)
    v = (3 * K - 2 * mu) / (2 * (3 * K + mu))
    return {"E": E, "v": v}, k_source


def febio_sheet_mixture(sheet_fit, matrix_fit, thickness_mm, warp_axis="x", weft_axis="y",
                        stiffness_scale=1.0, bulk_to_shear=100.0):
    """FEBio "solid mixture" for a sheet cast into a matrix, as
    (type, [(constituent type, {param: value}), ...], info)."""
    if sheet_fit.model not in SHEET_MODELS:
        raise UnsupportedMaterial(f"no FEBio sheet mapping for {sheet_fit.model!r}")
    if not thickness_mm or thickness_mm <= 0:
        raise ValueError("sheet layer needs a positive modelled thickness")
    p = sheet_fit.params
    lam0, beta = float(p["toe_stretch"]), float(p["toe_exponent"])
    c = fiber_pow_linear_c(lam0, beta)
    matrix, k_source = coupled_neo_hookean(matrix_fit, bulk_to_shear)
    parts = [("neo-Hookean", matrix)]
    info = {"matrix_E_MPa": matrix["E"], "matrix_v": matrix["v"], "matrix_k_source": k_source,
            "thickness_mm": thickness_mm, "stiffness_scale": stiffness_scale, "c": c}
    for name, S, axis in zip(("warp", "weft"), sheet_fit.warp_weft_N_per_mm(),
                             (warp_axis, weft_axis)):
        S = S * stiffness_scale
        E = S / (thickness_mm * c)
        parts.append(("fiber-pow-linear", {"E": E, "beta": beta, "lam0": lam0,
                                           "fiber": ",".join(f"{v:g}" for v in AXIS_VEC[axis])}))
        info[f"{name}_N_per_mm"] = S
        info[f"{name}_fiber_E_MPa"] = E
        info[f"{name}_axis"] = axis
    return "solid mixture", parts, info
