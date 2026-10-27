"""Map the curvature perception reports for a hand finger to a robot target.

Two modes, chosen in the finger TOML:
- "curvature" (default): kappa* = kappa_h. The robot finger copies the human
  finger's curvature.
- "angle": theta* = kappa_h * L_h, then kappa* = theta* / L_robot. The robot
  finger copies the bend angle. The robot finger (112 mm) is longer than a
  human finger (about 80-90 mm), so a direct curvature match bends it by
  roughly 25-40% more angle than this mode does.

The result is a curvature in 1/m. Clipping to what the robot can reach is the
shape map's job, not this module's.
"""

from dataclasses import dataclass

MODES = ("curvature", "angle")


@dataclass
class Mapping:
    mode: str = "curvature"
    L_h_mm: float = 85.0
    L_robot_mm: float = 112.0

    def __post_init__(self):
        if self.mode not in MODES:
            raise ValueError(f"mapping mode must be one of {MODES}, got {self.mode!r}")

    @classmethod
    def from_cfg(cls, cfg):
        return cls(cfg.mode, cfg.L_h_mm, cfg.L_robot_mm)

    def target_curvature(self, kappa_h):
        """Robot target curvature (1/m) for a human finger curvature (1/m)."""
        kappa_h = max(0.0, float(kappa_h))
        if self.mode == "curvature":
            return kappa_h
        theta = kappa_h * self.L_h_mm / 1000.0          # rad
        return theta / (self.L_robot_mm / 1000.0)
