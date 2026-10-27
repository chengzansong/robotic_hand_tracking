"""Closed-loop runner and step-response metrics, for tests and the demo.

`run_closed_loop` works with any lockstep Plant (see core/plant.py). Ground
truth columns (`P_true`, `kappa_true`) are filled only when the plant exposes
`true_pressure` and `true_curvature`, as LumpedPlant does.
"""

import numpy as np


def run_closed_loop(plant, finger, schedule, t_end, rate_hz=40.0):
    """Run `finger` on `plant`.

    schedule: list of (t_s, kappa_h) steps; the last one at or before t holds.
    Returns a dict of numpy arrays: t, P_true, P_hat, P_meas, P_target, kappa_true,
    kappa_target, duty, fault.
    """
    dt = 1.0 / rate_hz
    steps = sorted(schedule)
    rows = []
    while plant.now() < t_end - 1e-9:
        t = plant.now()
        kappa_h = [k for ts, k in steps if ts <= t + 1e-9][-1] if steps[0][0] <= t else 0.0
        samples = plant.read_pressure()
        duty = finger.tick(t, samples, kappa_h)
        plant.set_duty(finger.name, duty)
        est = finger.estimator
        raw = samples[-1].kPa if samples else np.nan
        rows.append((t, getattr(plant, "true_pressure", np.nan), est.P_hat,
                     est.cal.apply(raw) if samples else np.nan, finger.P_target,
                     getattr(plant, "true_curvature", np.nan), finger.kappa_target,
                     duty, int(est.fault)))
        plant.advance(dt)
    cols = ["t", "P_true", "P_hat", "P_meas", "P_target", "kappa_true", "kappa_target",
            "duty", "fault"]
    return dict(zip(cols, np.array(rows).T))


def step_metrics(t, y, t_step, t_next, target, band):
    """Settle time, overshoot and steady error of one step of y toward `target`.

    t, y: time and the measured quantity. The step starts at t_step and is
    judged until t_next. Returns a dict:
    - settle_s: time from t_step until y stays within +-band of target, or
      None if it never does before t_next.
    - overshoot: peak excursion beyond target in the direction of travel, in
      y's units (0 if none).
    - steady_err: mean(y) - target over the last 20% of the step.
    """
    m = (t >= t_step) & (t < t_next)
    tt, yy = t[m], y[m]
    y0 = yy[0]
    outside = np.abs(yy - target) > band
    if outside[-1]:
        settle = None
    elif not outside.any():
        settle = 0.0
    else:
        settle = float(tt[np.nonzero(outside)[0][-1] + 1] - t_step)
    direction = np.sign(target - y0) if target != y0 else 0.0
    over = float(max(0.0, np.max(direction * (yy - target)))) if direction else 0.0
    tail = tt >= t_next - 0.2 * (t_next - t_step)
    return {"settle_s": settle, "overshoot": over, "steady_err": float(yy[tail].mean() - target)}
