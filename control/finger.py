"""FingerController: one finger's whole control stack behind a single `tick`.

    perception curvature --mapping--> kappa* --shape map--> P*
    raw samples --estimator--> P_hat, dP_hat
    (t, P_hat, dP_hat, P*) --controller--> u --safety--> duty

It talks to no plant. The caller feeds it samples and sends the returned duty
wherever it goes (LumpedPlant, the future simulator, or the ESP32 link).
"""

from .controllers import Safety, build_controller
from .core.config import load_config
from .estimators.pressure import PressureEstimator
from .mapping import Mapping
from .shape.pressure_map import ShapeMap


class FingerController:
    def __init__(self, cfg, shape_map=None):
        self.cfg = cfg
        self.name = cfg.finger.name
        self.map = shape_map or ShapeMap.load(cfg.shape_map_path(), cfg.safety.p_cap_kPa)
        self.mapping = Mapping.from_cfg(cfg.mapping)
        self.estimator = PressureEstimator(cfg.estimator)
        self.controller = build_controller(cfg, self.map)
        self.safety = Safety(cfg.safety)
        # last tick, for logging
        self.kappa_target = self.P_target = self.duty = 0.0

    @classmethod
    def from_toml(cls, path):
        return cls(load_config(path))

    def tick(self, t, samples, kappa_h):
        """One control step. `kappa_h` is the human finger curvature in 1/m."""
        est = self.estimator
        est.update(samples)
        if samples:
            self.safety.feed(t)
        self.kappa_target = self.mapping.target_curvature(kappa_h)
        self.P_target = self.map.pressure_for(self.kappa_target)
        u = 0.0
        if est.ready:
            u = self.controller.step(t, est.P_hat, est.dP_hat, self.P_target)
        self.duty = self.safety.apply(t, u, est.P_hat, est.fault, est.ready)
        return self.duty
