"""Safety layer: always on, applied after whichever controller is selected.

- Pressure cap: at P_hat >= cap the duty is 0 until P_hat falls `cap_release`
  below the cap (a latch, so it cannot chatter at the cap). About 5.5 kPa for
  Ecoflex 00-30, since the finger curls into itself past about 6 kPa.
- Maximum duty: the duty is clamped to `max_duty`.
- Watchdog: if no fresh sensor sample has been fed for `watchdog_s`, u = 0.
  (The plant has its own watchdog on commands, as the firmware will.)
- Faults, or no startup zero yet: u = 0.
"""

from ..core.types import Fault


class Safety:
    def __init__(self, cfg):
        self.cfg = cfg
        self.reset()

    def reset(self):
        self.capped = False
        self.t_fed = None
        self.reason = ""

    def feed(self, t):
        """Call whenever a fresh sensor sample arrives."""
        self.t_fed = t

    def apply(self, t, u, P_hat, fault=Fault.NONE, ready=True):
        """Return the duty that is allowed. `reason` says why it was cut, if it was."""
        c = self.cfg
        self.reason = ""
        if self.t_fed is None or t - self.t_fed > c.watchdog_s:
            self.reason = "watchdog"
            return 0.0
        if fault:
            self.reason = f"fault {Fault(fault).name}"
            return 0.0
        if not ready:
            self.reason = "not zeroed"
            return 0.0
        if P_hat >= c.p_cap_kPa:
            self.capped = True
        elif P_hat < c.p_cap_kPa - c.cap_release_kPa:
            self.capped = False
        if self.capped:
            self.reason = "pressure cap"
            return 0.0
        return min(max(u, 0.0), c.max_duty)
