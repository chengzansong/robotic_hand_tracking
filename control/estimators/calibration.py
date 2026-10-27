"""Sensor calibration: P = gain * raw + offset, with the offset zeroed at startup.

`gain` comes from a bench calibration against a reference gauge (placeholder
1.0 until bench B6). The offset is measured at startup with the finger vented
(pump off, bleed open), by averaging raw readings for `zero_s` seconds.
"""


class Calibration:
    def __init__(self, gain=1.0, offset_kPa=0.0, zero_s=1.0):
        self.gain, self.offset_kPa, self.zero_s = gain, offset_kPa, zero_s
        self.zeroed = zero_s <= 0
        self._acc, self._n, self._t0 = 0.0, 0, None

    def apply(self, raw):
        return self.gain * raw + self.offset_kPa

    def zero_step(self, t, raw):
        """Feed one raw reading while vented. Returns True once the zero is set."""
        if self.zeroed:
            return True
        if self._t0 is None:
            self._t0 = t
        self._acc += raw
        self._n += 1
        if t - self._t0 >= self.zero_s:
            self.offset_kPa = -self.gain * self._acc / self._n
            self.zeroed = True
        return self.zeroed

    def rezero(self):
        """Start the zero routine again (finger must be vented)."""
        self.zeroed = self.zero_s <= 0
        self._acc, self._n, self._t0 = 0.0, 0, None
