"""LumpedPlant: one finger as a single-node ODE, for tests and the demo.

    duty u -> motor lag -> pump Q_p(u, P) -> [finger + short tube] -> bleed k*P^n
                                                   |
                                    shape map: curvature kappa(P)
                                                   |
                          sensor: gain, offset, noise, quantisation, 10/40 Hz, latency

    C(P) dP/dt = Q_p - Q_bleed        C(P) = V/P_atm + (1 + P/P_atm) dV/dP

C(P) comes from the FEM cavity volume V(P) (control/shape), plus the small
tube volume as extra gas. The tube's flow resistance is dropped: it is under
0.1 kPa for a hard tube shorter than 1 ft. The shape map is quasi-static, so
the finger has no bending lag. Integration is explicit Euler on a small fixed
step (1 ms), clamped at P >= 0.

The pump numbers are placeholders until the bench tests B2-B4 (see
control/README.md). Seeding makes every sensor draw repeatable.
"""

import numpy as np

from ..core import pneumatics as pn
from ..core.types import Sample


class LumpedPlant:
    def __init__(self, shape_map, pump, bleed, plant, sensor, finger="Index"):
        self.map = shape_map
        self.pump, self.bleed, self.cfg, self.sensor_cfg = pump, bleed, plant, sensor
        self.finger = finger
        self.bleed_k = bleed.k_ml_s if bleed.k_ml_s is not None else pn.orifice_k(
            bleed.diameter_mm, bleed.cd)
        self.P = 0.0                # true gauge pressure, kPa
        self.s = 0.0                # motor drive state, 0-1
        self.t = 0.0
        self.duty = 0.0
        self._t_cmd = 0.0
        self.stuck = False          # test hook: sensor repeats its last value
        rng = np.random.default_rng(plant.seed)
        self._rng = rng
        self.sensor_gain = 1.0 + rng.normal(0.0, sensor.gain_error_sigma)
        self.sensor_offset_kPa = rng.normal(0.0, sensor.offset_sigma_kPa)
        self._next_sample_t = 0.0
        self._seq = 0
        self._last_raw = 0.0
        self._pending = []          # (available_time, Sample)

    @classmethod
    def from_config(cls, cfg, shape_map, bleed_mm=None):
        """Build from a Config. `bleed_mm` overrides the orifice diameter."""
        bleed = cfg.bleed
        if bleed_mm is not None:
            bleed = type(bleed)(**{**vars(bleed), "diameter_mm": bleed_mm, "k_ml_s": None})
        return cls(shape_map, cfg.pump, bleed, cfg.plant, cfg.sensor, cfg.finger.name)

    # --- Plant protocol ---------------------------------------------------

    def set_duty(self, finger, duty):
        if finger != self.finger:
            raise KeyError(f"this plant drives {self.finger!r}, not {finger!r}")
        self.duty = float(np.clip(duty, 0.0, 1.0))
        self._t_cmd = self.t

    def read_pressure(self):
        ready = [s for ta, s in self._pending if ta <= self.t + 1e-12]
        self._pending = [(ta, s) for ta, s in self._pending if ta > self.t + 1e-12]
        return ready

    def now(self):
        return self.t

    def advance(self, dt_s):
        dt = self.cfg.dt_s
        n = max(1, round(dt_s / dt))
        for _ in range(n):
            self._step(dt_s / n)

    # --- ground truth, for scoring only (control must not read these) -------

    @property
    def true_pressure(self):
        return self.P

    @property
    def true_curvature(self):
        return float(self.map.kappa(self.P))

    # --- internals ----------------------------------------------------------

    def _step(self, dt):
        u = self.duty
        wd = self.cfg.watchdog_s
        if wd > 0 and self.t - self._t_cmd > wd:     # firmware-style command watchdog
            u = 0.0
        h = pn.duty_to_drive(u, self.pump.u_min)
        self.s += (h - self.s) * min(1.0, dt / max(self.pump.tau_motor_s, 1e-9))
        q_in = pn.pump_flow(self.s, self.P, self.pump.q0_ml_s, self.pump.p_stall_kPa)
        q_out = pn.bleed_flow(self.bleed_k, self.bleed.n, self.P)
        C = self.map.compliance_ml_per_kPa(self.P, self.cfg.v_tube_ml)
        self.P = max(0.0, self.P + dt * (q_in - q_out) / C)
        self.t += dt
        if self.t + 1e-12 >= self._next_sample_t:
            self._sample()
            self._next_sample_t += 1.0 / self.sensor_cfg.rate_hz

    def _sample(self):
        sc = self.sensor_cfg
        if self.stuck:
            raw = self._last_raw
        else:
            raw = self.sensor_gain * self.P + self.sensor_offset_kPa
            raw += self._rng.normal(0.0, sc.noise_sigma_kPa)
            raw = float(np.clip(raw, -sc.range_kPa, sc.range_kPa))
            raw = round(raw / sc.lsb_kPa) * sc.lsb_kPa
        self._last_raw = raw
        t_s = self._next_sample_t
        self._pending.append((t_s + sc.latency_s, Sample(t_s, self.finger, raw, self._seq)))
        self._seq += 1
