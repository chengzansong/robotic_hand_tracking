"""Closed-loop behaviour on LumpedPlant: pi_ff, bang_bang, open loop and faults."""

import numpy as np

from control.core.types import Fault, Sample
from control.estimators.pressure import PressureEstimator
from control.finger import FingerController
from control.runner import run_closed_loop, step_metrics

T0 = 1.5     # the sensor zero takes the first second or so


def run(finger, plant, kappas, hold_s=8.0, rate=40.0):
    sched = [(T0 + i * hold_s, k) for i, k in enumerate(kappas)]
    return run_closed_loop(plant, finger, sched, T0 + len(kappas) * hold_s, rate), sched


def test_pi_ff_settles_within_band_and_below_cap(make_loop, cfg):
    finger, plant = make_loop()
    log, sched = run(finger, plant, [17.0, 38.0, 27.0, 0.0])
    for i in (0, 1, 2):
        P = finger.map.pressure_for(sched[i][1])
        m = step_metrics(log["t"], log["P_true"], sched[i][0], sched[i][0] + 8.0, P, 0.3)
        assert m["settle_s"] is not None and m["settle_s"] < 3.0
        assert abs(m["steady_err"]) < 0.15
    assert log["P_true"].max() < cfg.safety.p_cap_kPa
    assert log["duty"].max() <= cfg.safety.max_duty + 1e-9


def test_pi_ff_never_exceeds_cap_at_the_largest_target(make_loop, cfg):
    finger, plant = make_loop()
    log, _ = run(finger, plant, [200.0], hold_s=10.0)      # far above reachable
    assert log["P_target"].max() <= finger.map.p_max + 1e-9
    assert log["P_true"].max() < cfg.safety.p_cap_kPa


def test_pi_ff_holds_against_a_bigger_bleed(make_loop):
    """The plant bleed is 0.4 mm but the controller believes 0.3 mm: it must still hold."""
    finger, plant = make_loop(bleed_mm=0.4)
    log, _ = run(finger, plant, [27.0], hold_s=15.0)
    P = finger.map.pressure_for(27.0)
    late = log["t"] > T0 + 8.0
    assert np.abs(log["P_true"][late] - P).max() < 0.35
    assert log["duty"][late].mean() > 0.4          # it is working against the leak


def test_pi_ff_learns_feedforward(make_loop):
    finger, plant = make_loop(bleed_mm=0.4)
    before = finger.controller.ff(3.0)
    run(finger, plant, [27.0], hold_s=20.0)
    assert finger.controller.ff(3.0) > before + 0.02


def test_bang_bang_stays_in_its_band(make_loop, cfg):
    finger, plant = make_loop(controller="bang_bang", ideal_sensor=True)
    log, _ = run(finger, plant, [27.0], hold_s=12.0)
    P, h = finger.map.pressure_for(27.0), cfg.controller.params["h_kPa"]
    late = log["t"] > T0 + 6.0
    # the band applies to what the controller sees; sensor latency, filter lag and one tick add up to about 0.3 kPa
    assert log["P_true"][late].min() > P - h - 0.3
    assert log["P_true"][late].max() < P + h + 0.3


def test_open_loop_volume_lands_near_target(make_loop):
    finger, plant = make_loop(controller="open_loop_volume", ideal_sensor=True)
    log, _ = run(finger, plant, [27.0], hold_s=10.0)
    assert abs(log["P_true"][-1] - finger.map.pressure_for(27.0)) < 0.8


def test_stuck_sensor_cuts_the_pump(make_loop):
    finger, plant = make_loop()
    log, _ = run(finger, plant, [27.0], hold_s=6.0)
    assert log["duty"][-1] > 0
    plant.stuck = True
    log = run_closed_loop(plant, finger, [(0.0, 27.0)], plant.now() + 2.0)
    assert finger.estimator.fault & Fault.STUCK
    assert log["duty"][-1] == 0.0
    assert plant.duty == 0.0


def test_out_of_range_sample_gives_zero_duty(cfg):
    f = FingerController(cfg)
    t = 0.0
    for i in range(80):                      # zero with clean samples
        t = i * 0.025
        f.tick(t, [Sample(t, "Index", 0.01 * (i % 3), i)], 20.0)
    assert f.estimator.ready and f.duty > 0
    t += 0.025
    u = f.tick(t, [Sample(t, "Index", 30.0, 99)], 20.0)
    assert f.estimator.fault & Fault.OUT_OF_RANGE
    assert u == 0.0
    t += 0.025
    assert f.tick(t, [Sample(t, "Index", 0.01, 100)], 20.0) == 0.0     # the fault latches


def test_impossible_slope_flagged(cfg):
    cfg.estimator.zero_s = 0.0
    est = PressureEstimator(cfg.estimator)
    for i, v in enumerate([0.0, 0.0, 0.0, 3.0, 3.0, 3.0]):      # 3 kPa in 25 ms, held
        est.update([Sample(i * 0.025, "Index", v + 0.01 * i, i)])
    assert est.fault & Fault.IMPOSSIBLE_SLOPE


def test_single_spike_is_rejected_by_the_median(cfg):
    cfg.estimator.zero_s = 0.0
    est = PressureEstimator(cfg.estimator)
    for i, v in enumerate([1.00, 1.01, 1.02, 4.0, 1.03, 1.04, 1.05]):
        est.update([Sample(i * 0.025, "Index", v, i)])
    assert est.fault == Fault.NONE
    assert abs(est.P_hat - 1.04) < 0.1
