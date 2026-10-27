"""The safety layer and the watchdogs."""

from control.controllers import Safety
from control.core.types import Fault


def test_fault_forces_zero(cfg):
    s = Safety(cfg.safety)
    s.feed(0.0)
    assert s.apply(0.01, 0.5, 2.0) == 0.5
    assert s.apply(0.02, 0.5, 2.0, fault=Fault.STUCK) == 0.0
    assert "STUCK" in s.reason


def test_watchdog_cuts_on_stale_samples_and_recovers(cfg):
    s = Safety(cfg.safety)
    assert s.apply(0.0, 0.5, 2.0) == 0.0               # nothing fed yet
    s.feed(0.0)
    assert s.apply(0.1, 0.5, 2.0) == 0.5
    assert s.apply(0.3, 0.5, 2.0) == 0.0               # 0.3 s without a sample
    assert s.reason == "watchdog"
    s.feed(0.31)
    assert s.apply(0.32, 0.5, 2.0) == 0.5


def test_cap_latches_and_max_duty_clamps(cfg):
    s = Safety(cfg.safety)
    s.feed(0.0)
    cap = cfg.safety.p_cap_kPa
    assert s.apply(0.0, 1.0, 1.0) == cfg.safety.max_duty
    assert s.apply(0.0, 0.5, cap + 0.1) == 0.0
    assert s.apply(0.0, 0.5, cap - 0.1) == 0.0          # still latched
    assert s.apply(0.0, 0.5, cap - 0.5) == 0.5          # released


def test_plant_watchdog_stops_the_pump_without_commands(make_loop):
    finger, plant = make_loop()
    plant.set_duty("Index", 0.9)
    plant.advance(0.15)
    assert plant.true_pressure > 0.05
    plant.advance(1.0)                                  # no new command: pump cut at 0.2 s
    p_a = plant.true_pressure
    plant.advance(0.5)
    assert plant.s < 0.01 and plant.true_pressure < p_a    # draining through the bleed
