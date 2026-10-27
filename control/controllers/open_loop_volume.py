"""open_loop_volume: run the pump for as long as the volume needs, no sensor.

When the target rises, the pump runs at u_run for t = dV / Q_p, where dV is
the extra standard volume the FEM cavity needs, (P_atm + P) V(P) / P_atm
between the old and new target, and Q_p is the believed pump flow at the
mid pressure. It then switches to the model hold duty (the duty that balances
the believed bleed), or off if hold = false. When the target falls, it stops.
Use it for sanity checks or when the sensor is out. It ignores P_hat, so any
error in the pump or leak numbers goes straight into the pressure.
"""

import numpy as np

from ..core import pneumatics as pn
from ..shape.pressure_map import P_ATM_KPA
from .pi_ff import model_ff_table


class OpenLoopVolume:
    def __init__(self, shape_map, model, p_max_kPa, u_run=0.8, hold=True, **_unused):
        self.map, self.m, self.u_run, self.hold = shape_map, model, u_run, hold
        self.p_grid, self.u_hold = model_ff_table(model, p_max_kPa)
        self.reset()

    def reset(self):
        self.p_prev = 0.0          # pressure the finger is believed to be at
        self.t_end = -1.0
        self.t_last_target = None

    def _std_volume(self, P):
        return (P_ATM_KPA + P) * float(self.map.volume_ml(P)) / P_ATM_KPA

    def step(self, t, P_hat, dP_hat, P_target):
        if self.t_last_target is None or abs(P_target - self.t_last_target) > 1e-6:
            if P_target > self.p_prev:
                dv = self._std_volume(P_target) - self._std_volume(self.p_prev)
                h = pn.duty_to_drive(self.u_run, self.m.u_min)
                q = pn.pump_flow(h, 0.5 * (P_target + self.p_prev), self.m.q0_ml_s,
                                 self.m.p_stall_kPa)
                self.t_end = t + dv / q if q > 0 else t
            else:
                self.t_end = t
            self.p_prev = P_target
            self.t_last_target = P_target
        if P_target <= 0.05:
            return 0.0
        if t < self.t_end:
            return self.u_run
        return float(np.interp(P_target, self.p_grid, self.u_hold)) if self.hold else 0.0
