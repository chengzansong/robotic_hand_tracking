"""Refit Marechal et al. 2021 raw uniaxial data in the canonical (Abaqus) conventions.

Ogden  P = sum 2 mu_i/alpha_i (lam^(alpha_i-1) - lam^(-alpha_i/2-1)),  mu_i >= 0 enforced (each term of W = 2 mu_i/alpha_i^2 (...) is then non-negative and stable)
Yeoh   P = 2 (lam - lam^-2) (C10 + 2 C20 (I1-3) + 3 C30 (I1-3)^2),       C10 >= 0, C30 >= 0
Nominal (engineering) stress vs stretch, incompressible uniaxial tension.

Usage: python -I fit.py <Tensile-Tests-Data dir> [lam_max]
"""
import itertools
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear, nnls

DATA = Path(sys.argv[1])
LAM_MAX = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
PSI = 6.894757

FILES = {
    "ecoflex_00-10": ("Ecoflex 00-10.csv", 8),
    "ecoflex_00-30": ("Ecoflex 00-30.csv", 10),
    "ecoflex_00-50": ("Ecoflex 00-50.csv", 12),
    "dragon_skin_10": ("Dragon Skin 10 MEDIUM.csv", 22),
    "dragon_skin_20": ("Dragon Skin 20.csv", 49),
    "dragon_skin_30": ("Dragon Skin 30.csv", 86),
}


def load(name):
    rows, started = [], False
    for line in (DATA / name).read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("Time (s)"):
            started = True
            continue
        if started and line.strip():
            try:
                rows.append([float(x) for x in line.split(";")[:5]])
            except ValueError:
                pass
    a = np.array(rows)
    ee, se = a[:, 3], a[:, 4] * 1000.0  # kPa
    keep = (ee > 0) & np.isfinite(se)
    ee, se = ee[keep], se[keep]
    o = np.argsort(ee)
    return 1 + ee[o], se[o]


def ogden_basis(lam, alpha):
    return 2 / alpha * (lam ** (alpha - 1) - lam ** (-alpha / 2 - 1))


def ogden_P(mu, al, lam):
    return sum(m * ogden_basis(lam, a) for m, a in zip(mu, al))


def ogden_eqb(mu, al, lam):
    # equibiaxial nominal stress, incompressible: P = sum 2mu/a (lam^(a-1) - lam^(-2a-1))
    return sum(2 * m / a * (lam ** (a - 1) - lam ** (-2 * a - 1)) for m, a in zip(mu, al))


def yeoh_P(c, lam):
    I1m3 = lam ** 2 + 2 / lam - 3
    return 2 * (lam - lam ** -2) * (c[0] + 2 * c[1] * I1m3 + 3 * c[2] * I1m3 ** 2)


def yeoh_eqb(c, lam):
    I1m3 = 2 * lam ** 2 + lam ** -4 - 3
    return 2 * (lam - lam ** -5) * (c[0] + 2 * c[1] * I1m3 + 3 * c[2] * I1m3 ** 2)


def fit_ogden(lam, P, w, n=3):
    grid = [a for a in np.arange(-10, 12.01, 0.5) if abs(a) > 1e-9]
    best = None
    for als in itertools.combinations(grid, n):
        A = np.column_stack([ogden_basis(lam, a) for a in als]) * w[:, None]
        x, r = nnls(A, P * w)
        if best is None or r < best[0]:
            best = (r, als, x)
    # local refinement on a finer grid
    r0, als0, _ = best
    for step in (0.1, 0.02):
        improved = True
        while improved:
            improved = False
            for i in range(n):
                for d in (-step, step):
                    als = list(best[1])
                    als[i] = round(als[i] + d, 4)
                    if abs(als[i]) < 1e-6:
                        continue
                    A = np.column_stack([ogden_basis(lam, a) for a in als]) * w[:, None]
                    x, r = nnls(A, P * w)
                    if r < best[0] - 1e-12:
                        best = (r, tuple(als), x)
                        improved = True
    r, als, x = best
    mu = list(x)
    # drop terms with zero weight
    pairs = [(m, a) for m, a in zip(mu, als) if abs(m) > 1e-9]
    return [p[0] for p in pairs], [p[1] for p in pairs]


def fit_yeoh(lam, P, w):
    I1m3 = lam ** 2 + 2 / lam - 3
    f = 2 * (lam - lam ** -2)
    A = np.column_stack([f, f * 2 * I1m3, f * 3 * I1m3 ** 2]) * w[:, None]
    res = lsq_linear(A, P * w, bounds=([0, -np.inf, 0], [np.inf, np.inf, np.inf]))
    return list(res.x)


def check_stable(fn, lo=0.3, hi=6.0):
    lam = np.linspace(lo, hi, 2000)
    P = fn(lam)
    return bool(np.all(np.diff(P) > 0))


for mid, (fn, ds_psi) in FILES.items():
    lam_all, P_all = load(fn)
    m = lam_all <= LAM_MAX
    lam_d, P_d = lam_all[m], P_all[m]
    # resample uniformly in stretch so dense regions do not dominate
    lam = np.linspace(1.02, min(LAM_MAX, lam_d.max()), 150)
    P = np.interp(lam, lam_d, P_d)
    w = 1.0 / (P + 0.1 * P.max())  # blend of relative and absolute error
    mu, al = fit_ogden(lam, P, w)
    c = fit_yeoh(lam, P, w)
    ds = ds_psi * PSI
    data2 = float(np.interp(2.0, lam_d, P_d))
    print(f"\n== {mid}  fit range lam 1..{lam.max():.2f}  datasheet P100 = {ds:.1f} kPa  data P(2) = {data2:.1f} kPa")
    print(f"  OGDEN mu_kPa = {[float(f'{v:.5g}') for v in mu]}  alpha = {[float(f'{v:.4g}') for v in al]}")
    print(f"        mu0 = {sum(mu):.2f} kPa  P(2) = {ogden_P(mu, al, 2.0):.1f}  ratio fit/ds = {ogden_P(mu, al, 2.0)/ds:.2f}  "
          f"eqbiax P(2) = {ogden_eqb(mu, al, 2.0):.1f}  stable(uniax 0.3-6) = {check_stable(lambda l: ogden_P(mu, al, l))}  stable(eqbiax 1-4) = {check_stable(lambda l: ogden_eqb(mu, al, l), 1.0, 4.0)}")
    print(f"  YEOH  c10,c20,c30 kPa = {[float(f'{v:.5g}') for v in c]}")
    print(f"        mu0 = {2*c[0]:.2f} kPa  P(2) = {yeoh_P(c, 2.0):.1f}  ratio fit/ds = {yeoh_P(c, 2.0)/ds:.2f}  "
          f"eqbiax P(2) = {yeoh_eqb(c, 2.0):.1f}  stable(uniax 0.3-6) = {check_stable(lambda l: yeoh_P(c, l))}  "
          f"stable(eqbiax 1-4) = {check_stable(lambda l: yeoh_eqb(c, l), 1.0, 4.0)}")
    for name, f in (("ogden", lambda l: ogden_P(mu, al, l)), ("yeoh", lambda l: yeoh_P(c, l))):
        errs = []
        for lq in (1.1, 1.25, 1.5, 2.0, 3.0, 4.0):
            if lq <= lam.max():
                pd_ = float(np.interp(lq, lam_d, P_d))
                errs.append(f"{lq}:{100*(f(lq)-pd_)/pd_:+.0f}%")
        rms = np.sqrt(np.mean((f(lam) - P) ** 2))
        print(f"        {name} rel err vs data at lam " + " ".join(errs) + f"   RMS {rms:.2f} kPa")
