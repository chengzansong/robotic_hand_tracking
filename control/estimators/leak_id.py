"""Leak identification: fit the bleed law Q = k * P^n from a pressure decay.

With the pump off, C(P) dP/dt = -k P^n. The fit integrates that ODE from a
fitted start pressure and matches it to the measured P(t) by least squares, so
it never differentiates noisy data. C(P) is the controller's belief, from the
FEM cavity volume. k is in ml/s per kPa^n.
"""

import numpy as np
from scipy.optimize import least_squares


def simulate_decay(t, P0, k, n, compliance, dt=0.002):
    """Integrate C(P) dP/dt = -k P^n on a fine grid and sample it at times t."""
    t = np.asarray(t, float) - t[0]
    tg = np.arange(0.0, t[-1] + dt, dt)
    P = np.empty_like(tg)
    p = max(P0, 0.0)
    for i in range(len(tg)):
        P[i] = p
        p = max(0.0, p - dt * k * p**n / compliance(p))
    return np.interp(t, tg, P)


def fit_leak(t, P, compliance, n=0.5, fit_n=False, p_floor=0.3):
    """Fit k (and optionally n) to a decay trace.

    t, P: sample times (s) and calibrated pressures (kPa) after the pump is off.
    compliance: callable P -> ml/kPa.
    Points below `p_floor` kPa are ignored: sensor noise dominates there.
    Returns {"k": ..., "n": ..., "rms_kPa": ...}.
    """
    t, P = np.asarray(t, float), np.asarray(P, float)
    keep = P > p_floor
    if keep.sum() < 5:
        raise ValueError("decay trace has too few points above the floor")
    # use only the leading run above the floor
    end = len(keep) if keep.all() else int(np.argmin(keep))
    t, P = t[:end], P[:end]

    def resid(x):
        nn = x[2] if fit_n else n
        return simulate_decay(t, x[1], np.exp(x[0]), nn, compliance) - P

    x0 = [0.0, float(P[:3].mean())] + ([n] if fit_n else [])
    lo = [-6.0, 0.0] + ([0.4] if fit_n else [])
    hi = [4.0, 15.0] + ([1.2] if fit_n else [])
    sol = least_squares(resid, x0, bounds=(lo, hi), diff_step=1e-3)
    nn = float(sol.x[2]) if fit_n else n
    return {"k": float(np.exp(sol.x[0])), "n": nn,
            "rms_kPa": float(np.sqrt(np.mean(sol.fun**2)))}


def run_decay_experiment(plant, finger, estimator, compliance, p_fill=3.5, u_fill=0.6,
                         t_fill_max=6.0, t_decay=3.0, t_skip=0.4, n=0.5, fit_n=False,
                         tick_s=0.025):
    """Pump to p_fill, switch off, record the decay and fit the leak.

    The first `t_skip` seconds after switch-off are dropped: the pump is still
    spinning down and the sensor lags, so the pressure is still rising there.

    Uses only the Plant protocol, and assumes the estimator is already zeroed.
    Returns the dict from fit_leak.
    """
    t0 = plant.now()
    while estimator.P_hat < p_fill and plant.now() - t0 < t_fill_max:
        plant.set_duty(finger, u_fill)
        plant.advance(tick_s)
        estimator.update(plant.read_pressure())
    ts, ps = [], []
    t1 = plant.now()
    while plant.now() - t1 < t_decay:
        plant.set_duty(finger, 0.0)
        plant.advance(tick_s)
        for s in plant.read_pressure():
            estimator.update([s])
            if plant.now() - t1 < t_skip:
                continue
            ts.append(s.t_s)
            ps.append(estimator.cal.apply(s.kPa))     # unfiltered, so filter lag cannot bias the fit
    return fit_leak(ts, ps, compliance, n=n, fit_n=fit_n)
