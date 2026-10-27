"""PressureEstimator: raw samples in, P_hat, dP_hat and fault flags out.

Chain per sample: calibration, median-of-3, then an alpha-beta (default) or
EMA filter. Faults latch until `clear_faults()`:
- OUT_OF_RANGE: calibrated reading outside [range_lo, range_hi].
- STUCK: `stuck_n` identical raw readings in a row (real noise is never that flat).
- IMPOSSIBLE_SLOPE: the median-filtered pressure moves faster than the pump
  could drive it, `max_slope_kPa_s`.
Until the startup zero is done, `ready` is False and the loop commands u = 0.
"""

from ..core.types import Fault
from .calibration import Calibration
from .filters import Median3, make_filter


class PressureEstimator:
    def __init__(self, cfg):
        self.cfg = cfg
        self.cal = Calibration(cfg.gain, cfg.offset_kPa, cfg.zero_s)
        self.median = Median3()
        self.filt = make_filter(cfg.filter, cfg)
        self.reset()

    def reset(self):
        """Clear filter state and faults. The calibration is kept."""
        self.P_hat = 0.0
        self.dP_hat = 0.0
        self.t_last = None
        self.fault = Fault.NONE
        self._prev_med = self._prev_t = None
        self._last_raw, self._same = None, 0
        self.median.reset()
        self.filt.reset()

    @property
    def ready(self):
        return self.cal.zeroed

    def clear_faults(self):
        self.fault = Fault.NONE
        self._same = 0

    def update(self, samples):
        """Process new raw Samples (oldest first). Returns the latest P_hat."""
        for s in samples:
            self._one(s)
        return self.P_hat

    def _one(self, s):
        c = self.cfg
        if not self.cal.zeroed:
            self.cal.zero_step(s.t_s, s.kPa)
            self._last_raw = s.kPa
            return
        if self._last_raw is not None and s.kPa == self._last_raw:
            self._same += 1
            if self._same >= c.stuck_n:
                self.fault |= Fault.STUCK
        else:
            self._same = 0
        self._last_raw = s.kPa
        p = self.cal.apply(s.kPa)
        if not (c.range_lo_kPa <= p <= c.range_hi_kPa):
            self.fault |= Fault.OUT_OF_RANGE
        med = self.median.update(p)
        if self._prev_med is not None:
            slope = abs(med - self._prev_med) / max(s.t_s - self._prev_t, 1e-6)
            if slope > c.max_slope_kPa_s:
                self.fault |= Fault.IMPOSSIBLE_SLOPE
        self._prev_med, self._prev_t = med, s.t_s
        self.P_hat, self.dP_hat = self.filt.update(s.t_s, med)
        self.t_last = s.t_s
