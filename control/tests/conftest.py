"""Shared fixtures: the default config and a loop builder."""

from pathlib import Path

import pytest

from control.core.config import load_config
from control.finger import FingerController
from control.plants.lumped import LumpedPlant

CONFIG = Path(__file__).resolve().parent.parent / "configs" / "finger_eco30.toml"


@pytest.fixture
def cfg():
    return load_config(CONFIG)


@pytest.fixture
def make_loop(cfg):
    """make_loop(bleed_mm, ideal_sensor, controller) -> (finger, plant)."""

    def build(bleed_mm=None, ideal_sensor=False, controller=None):
        if controller:
            cfg.controller.type = controller
        if ideal_sensor:
            cfg.sensor.noise_sigma_kPa = 0.02    # not zero: a perfectly flat sensor trips the stuck fault
            cfg.sensor.gain_error_sigma = 0.0
            cfg.sensor.offset_sigma_kPa = 0.0
        finger = FingerController(cfg)
        plant = LumpedPlant.from_config(cfg, finger.map, bleed_mm=bleed_mm)
        return finger, plant

    return build
