"""Calibration, filters, config loading and leak identification."""

import numpy as np
import pytest

from control.core.config import config_from_dict
from control.core.types import Sample
from control.estimators.calibration import Calibration
from control.estimators.filters import EMA, AlphaBeta
from control.estimators.leak_id import fit_leak, run_decay_experiment, simulate_decay
from control.estimators.pressure import PressureEstimator
from control.shape.pressure_map import ShapeMap


def test_startup_zero_removes_the_offset(cfg):
    est = PressureEstimator(cfg.estimator)
    for i in range(60):
        est.update([Sample(i * 0.025, "Index", 0.37, i)])
    assert est.ready
    assert est.cal.apply(2.37) == pytest.approx(2.0)


def test_calibration_gain():
    assert Calibration(gain=1.1, zero_s=0.0).apply(2.0) == pytest.approx(2.2)


def test_alpha_beta_tracks_a_ramp():
    f = AlphaBeta(0.5, 0.15)
    for i in range(80):
        x, v = f.update(i * 0.025, 2.0 * i * 0.025)
    assert v == pytest.approx(2.0, abs=0.1)
    assert x == pytest.approx(2.0 * 79 * 0.025, abs=0.1)


def test_ema_converges():
    f = EMA(0.4)
    for i in range(50):
        x, v = f.update(i * 0.025, 3.0)
    assert x == pytest.approx(3.0, abs=1e-3) and abs(v) < 0.01


def test_config_rejects_unknown_keys():
    with pytest.raises(ValueError):
        config_from_dict({"pump": {"q0_ml": 5.0}})


def test_leak_fit_recovers_known_k_from_a_noisy_decay(cfg):
    m = ShapeMap.load(cfg.shape_map_path(), cfg.finger.p_cap_kPa)

    def comp(P):
        return m.compliance_ml_per_kPa(P, 1.0)

    t = np.arange(0, 3.0, 0.025)
    P = simulate_decay(t, 3.5, 2.0, 0.5, comp)
    P = P + np.random.default_rng(3).normal(0, 0.08, len(t))
    assert fit_leak(t, P, comp)["k"] == pytest.approx(2.0, rel=0.08)
    assert 0.4 <= fit_leak(t, P, comp, fit_n=True)["n"] <= 1.2


def test_leak_fit_through_the_plant(make_loop, cfg):
    finger, plant = make_loop(bleed_mm=0.4)
    while not finger.estimator.ready:
        plant.advance(0.025)
        finger.estimator.update(plant.read_pressure())

    def comp(P):
        return finger.map.compliance_ml_per_kPa(P, cfg.plant.v_tube_ml)

    fit = run_decay_experiment(plant, "Index", finger.estimator, comp)
    assert fit["k"] == pytest.approx(plant.bleed_k, rel=0.12)
