"""The Plant protocol: what control needs from a finger and its pump.

`LumpedPlant` (control/plants/lumped.py) implements it for tests and the demo.
The real-time simulator in simulation/stage2_hand/ will implement the same
three calls (as `SimIO`) and the future ESP32 bridge will too. Control code
never knows which one it is talking to, and never imports simulation/.
"""

from typing import Protocol

from .types import Sample


class Plant(Protocol):
    def set_duty(self, finger: str, duty: float) -> None:
        """Pump duty 0-1 for one finger. Latched (zero-order hold) until changed."""

    def read_pressure(self) -> list[Sample]:
        """New raw pressure samples since the last call, oldest first."""

    def now(self) -> float:
        """Plant time in seconds (simulated or wall clock)."""

    def advance(self, dt_s: float) -> None:
        """Lockstep plants only: run the plant forward by dt_s.

        A real-time plant (hardware, or a simulator paced to the wall clock)
        makes this a no-op or a sleep.
        """
