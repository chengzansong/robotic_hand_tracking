"""pi_ff (default): learned feedforward plus a PI correction.

    u = u_ff(P*) + Kp e + Ki I - Kd dP_hat,        e = P* - P_hat

u_ff(P*) is the duty at which pump flow equals bleed flow at P*. It starts
from the controller's model of the pump and bleed (`model_ff_table`), and
while the pressure sits at the target it learns: the integrator's contribution
is moved into the table a little at a time, so the next visit to that pressure
starts with the right duty. Anti-windup is conditional integration, which
stops accumulating while the output is saturated in the direction of the
error, and the integrator is also clamped. The output is clamped to [0, u_max].

The pump only adds air, so with e < 0 the duty falls to 0 and the bleed
brings the pressure down. `pi_ff_valve` would open a release valve instead,
and is a stub until a valve exists.
"""

import numpy as np

from ..core import pneumatics as pn


def model_ff_table(model, p_max_kPa, step_kPa=0.5):
    """Duty that balances the believed bleed, on a pressure grid 0..p_max.

    Returns (p_grid, u_ff). The duty at P = 0 is 0 (no leak to balance).
    """
    k = model.k_ml_s if model.k_ml_s is not None else pn.orifice_k(model.bleed_mm, model.cd)
    grid = np.arange(0.0, p_max_kPa + 1e-9, step_kPa)
    u = np.zeros_like(grid)
    for i, P in enumerate(grid):
        q_pump_full = pn.pump_flow(1.0, P, model.q0_ml_s, model.p_stall_kPa)
        if P > 0 and q_pump_full > 0:
            h = pn.bleed_flow(k, model.n, P) / q_pump_full
            u[i] = min(1.0, pn.drive_to_duty(h, model.u_min))
    return grid, u


class PIFF:
    def __init__(self, shape_map=None, model=None, p_max_kPa=5.0, kp=0.7, ki=0.5,
                 kd=0.05, i_max=0.4, u_max=1.0, ff_learn=True, ff_rate_per_s=0.3,
                 settle_band_kPa=0.15, settle_slope_kPa_s=0.4, settle_hold_s=0.5,
                 ff_table=None, **_unused):
        self.kp, self.ki, self.kd = kp, ki, kd
        self.i_max, self.u_max = i_max, u_max
        self.ff_learn, self.ff_rate = ff_learn, ff_rate_per_s
        self.band, self.slope, self.hold_s = settle_band_kPa, settle_slope_kPa_s, settle_hold_s
        if ff_table is not None:                     # [[p_kPa, u], ...] from the TOML
            self.p_grid, self.u_ff = (np.array(c, float) for c in zip(*ff_table))
        else:
            self.p_grid, self.u_ff = model_ff_table(model, p_max_kPa)
        self.reset()

    def reset(self):
        self.I = 0.0
        self.t_prev = None
        self.t_settled = 0.0
        self.last_u_ff = 0.0

    def ff(self, P_target):
        return float(np.interp(P_target, self.p_grid, self.u_ff))

    def learned_table(self):
        """The current u_ff table as (p_grid, u_ff) arrays, for saving."""
        return self.p_grid.copy(), self.u_ff.copy()

    def step(self, t, P_hat, dP_hat, P_target):
        dt = 0.0 if self.t_prev is None else t - self.t_prev
        self.t_prev = t
        if P_target <= 0.05:                         # rest: no pump, no windup
            self.I = 0.0
            self.t_settled = 0.0
            return 0.0
        e = P_target - P_hat
        u_ff = self.last_u_ff = self.ff(P_target)
        u_raw = u_ff + self.kp * e + self.ki * self.I - self.kd * dP_hat
        u = min(max(u_raw, 0.0), self.u_max)
        saturated = (u_raw > self.u_max and e > 0) or (u_raw < 0.0 and e < 0)
        if not saturated:
            self.I = float(np.clip(self.I + e * dt, -self.i_max / self.ki if self.ki else 0.0,
                                   self.i_max / self.ki if self.ki else 0.0))
        if self.ff_learn:
            self._learn(P_target, e, dP_hat, saturated, dt)
        return u

    def _learn(self, P_target, e, dP_hat, saturated, dt):
        """Move the integrator's duty into the u_ff table while settled."""
        settled = abs(e) < self.band and abs(dP_hat) < self.slope and not saturated
        self.t_settled = self.t_settled + dt if settled else 0.0
        if self.t_settled < self.hold_s or self.ki == 0:
            return
        move = min(1.0, self.ff_rate * dt) * self.ki * self.I     # duty to transfer
        j = int(np.clip(np.searchsorted(self.p_grid, P_target) - 1, 0, len(self.p_grid) - 2))
        w = (P_target - self.p_grid[j]) / (self.p_grid[j + 1] - self.p_grid[j])
        w = float(np.clip(w, 0.0, 1.0))
        # keep the interpolated u_ff at P_target shifted by exactly `move`
        norm = (1 - w) ** 2 + w ** 2
        self.u_ff[j] += move * (1 - w) / norm
        self.u_ff[j + 1] += move * w / norm
        self.u_ff = np.clip(self.u_ff, 0.0, 1.0)
        self.I -= move / self.ki


class PIFFValve:
    """Stub: pi_ff plus a release valve that opens when e < -band. No valve yet."""

    def __init__(self, **_unused):
        raise NotImplementedError("pi_ff_valve needs a release valve, which v1 does not have")
