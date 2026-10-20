"""Independent check of every material card through the repo loader.

Usage: python validate_cards.py <repo root> <Marechal Tensile-Tests-Data dir>
"""
import sys
from pathlib import Path

import numpy as np

repo, data_dir = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(repo))
from simulation.core.materials import list_materials, load_material  # noqa: E402

RAW = {
    "ecoflex_00-10": "Ecoflex 00-10.csv", "ecoflex_00-30": "Ecoflex 00-30.csv",
    "ecoflex_00-50": "Ecoflex 00-50.csv", "dragon_skin_10": "Dragon Skin 10 MEDIUM.csv",
    "dragon_skin_20": "Dragon Skin 20.csv", "dragon_skin_30": "Dragon Skin 30.csv",
}


def P_uniax(fit, lam):
    p, m = fit.params, fit.model
    if m == "neo_hookean":
        return p["mu_kPa"] * (lam - lam ** -2)
    if m == "mooney_rivlin":
        return 2 * (lam - lam ** -2) * (p["c10_kPa"] + p["c01_kPa"] / lam)
    if m == "yeoh":
        I = lam ** 2 + 2 / lam - 3
        return 2 * (lam - lam ** -2) * (p["c10_kPa"] + 2 * p["c20_kPa"] * I + 3 * p["c30_kPa"] * I ** 2)
    if m == "ogden":
        return sum(2 * u / a * (lam ** (a - 1) - lam ** (-a / 2 - 1)) for u, a in zip(p["mu_kPa"], p["alpha"]))


def P_eqb(fit, lam):
    p, m = fit.params, fit.model
    if m == "neo_hookean":
        return p["mu_kPa"] * (lam - lam ** -5)
    if m == "yeoh":
        I = 2 * lam ** 2 + lam ** -4 - 3
        return 2 * (lam - lam ** -5) * (p["c10_kPa"] + 2 * p["c20_kPa"] * I + 3 * p["c30_kPa"] * I ** 2)
    if m == "ogden":
        return sum(2 * u / a * (lam ** (a - 1) - lam ** (-2 * a - 1)) for u, a in zip(p["mu_kPa"], p["alpha"]))
    return float("nan")


def raw_P(name, lam):
    rows, started = [], False
    for line in (data_dir / name).read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("Time (s)"):
            started = True
            continue
        if started and line.strip():
            try:
                rows.append([float(x) for x in line.split(";")[:5]])
            except ValueError:
                pass
    a = np.array(rows)
    e, s = a[:, 3], a[:, 4] * 1000
    o = np.argsort(e)
    return float(np.interp(lam - 1, e[o], s[o]))


for mid in list_materials():
    mat = load_material(mid)
    ds = mat.datasheet.get("modulus_100_kPa")
    print(f"\n{mid}  density {mat.density_kg_m3}  datasheet P100 {ds}")
    pur = __import__("tomllib").loads((repo / "simulation/materials" / f"{mid}.toml").read_text(encoding="utf-8")).get("purchase", {})
    if "price_usd" in pur:
        calc = pur["price_usd"] / pur["kit_mass_kg"]
        print(f"  purchase {pur['kit_size']}: ${pur['price_usd']} / {pur['kit_mass_kg']} kg = {calc:.2f} $/kg (card {pur['price_per_kg_usd']})"
              + ("  MISMATCH" if abs(calc - pur["price_per_kg_usd"]) > 0.05 else ""))
    else:
        print(f"  purchase: no USD price ({pur.get('vendor')})")
    for i, f in enumerate(mat.fits):
        p2 = P_uniax(f, 2.0)
        ratio = p2 / ds if ds else float("nan")
        lam = np.linspace(0.3, 6.0, 3000)
        mono = bool(np.all(np.diff(P_uniax(f, lam)) > 0))
        lamb = np.linspace(1.0, 4.0, 1000)
        mono_b = bool(np.all(np.diff(P_eqb(f, lamb)) > 0)) if f.model != "mooney_rivlin" else None
        flag = "  FLAG >2x" if ds and (ratio > 2 or ratio < 0.5) else ""
        extra = ""
        if mid in RAW and "Marechal" in f.source:
            errs = []
            for lq in (1.1, 1.5, 2.0, 3.0):
                d = raw_P(RAW[mid], lq)
                errs.append(f"{lq}:{100*(P_uniax(f, lq)-d)/d:+.0f}%")
            extra = "  vs Marechal raw " + " ".join(errs)
        print(f"  [{i}] {f.model:11s} {f.confidence:6s} mu0 {f.initial_shear_modulus_kPa():7.2f}  P(2) {p2:7.1f}  ratio {ratio:5.2f}"
              f"  eqb(2) {P_eqb(f, 2.0):7.1f}  stable uniax {mono} eqb {mono_b}{flag}{extra}")
