"""Small types shared across control."""

from dataclasses import dataclass
from enum import IntFlag


@dataclass(frozen=True)
class Sample:
    """One pressure reading as the controller receives it (raw sensor units).

    `kPa` is the raw gauge reading, before calibration. `t_s` is the time the
    reading was taken, not the time it arrived, so latency shows up as
    `now() - t_s`.
    """

    t_s: float
    finger: str
    kPa: float
    seq: int


class Fault(IntFlag):
    """Sensor fault flags raised by the estimator. Any set flag forces u = 0."""

    NONE = 0
    OUT_OF_RANGE = 1
    STUCK = 2
    IMPOSSIBLE_SLOPE = 4
