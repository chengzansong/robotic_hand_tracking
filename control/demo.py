"""Closed-loop demo: one finger on LumpedPlant, stepping through curvature targets.

    python -m control.demo
    python -m control.demo --controller bang_bang --bleed-mm 0.3

Prints settle time, overshoot and steady error in pressure and curvature for
each step, and saves a plot to control/outputs/ (gitignored). `--bleed-mm`
changes the plant's bleed only; the controller's belief stays as in the config,
so a different size tests how well the controller copes with a wrong leak model.
"""

import argparse
from pathlib import Path

from .controllers import build_controller
from .core.config import CONTROL_DIR, load_config
from .estimators.leak_id import run_decay_experiment
from .finger import FingerController
from .plants.lumped import LumpedPlant
from .runner import run_closed_loop, step_metrics

DEFAULT_CONFIG = CONTROL_DIR / "configs" / "finger_eco30.toml"
T_ZERO_S = 1.5          # time at the start for the sensor zero
T_STEP_S = 6.0          # time per target
TARGETS = [0.0, 17.0, 38.0, 27.0, 0.0]    # curvature, 1/m
P_BAND_KPA = 0.2
K_BAND = 2.0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", default=str(DEFAULT_CONFIG))
    ap.add_argument("--controller", help="override [controller] type")
    ap.add_argument("--bleed-mm", type=float, help="plant bleed diameter (mm)")
    ap.add_argument("--targets", type=float, nargs="+", default=TARGETS,
                    help="curvature targets in 1/m, one per step")
    ap.add_argument("--identify-leak", action="store_true",
                    help="run a decay test first and give the controller the fitted leak")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--out", default=str(CONTROL_DIR / "outputs" / "demo.png"))
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    if a.controller:
        cfg.controller.type = a.controller
    finger = FingerController(cfg)
    plant = LumpedPlant.from_config(cfg, finger.map, bleed_mm=a.bleed_mm)
    print(f"controller={cfg.controller.type}  plant bleed k={plant.bleed_k:.2f} ml/s/kPa^n "
          f"(belief: {cfg.model.bleed_mm} mm)  sensor offset={plant.sensor_offset_kPa:+.2f} kPa "
          f"gain={plant.sensor_gain:.3f}")

    # zero the sensor first (finger vented, pump off)
    while not finger.estimator.ready:
        plant.advance(1 / cfg.loop.rate_hz)
        finger.estimator.update(plant.read_pressure())
    if a.identify_leak:
        def comp(P):
            return finger.map.compliance_ml_per_kPa(P, cfg.plant.v_tube_ml)
        fit = run_decay_experiment(plant, finger.name, finger.estimator, comp, n=cfg.bleed.n)
        print(f"identified leak: k={fit['k']:.2f} ml/s/kPa^n (rms {fit['rms_kPa']:.3f} kPa)")
        cfg.model.k_ml_s, cfg.model.n = fit["k"], fit["n"]
        finger.controller = build_controller(cfg, finger.map)
        while plant.true_pressure > 0.05:        # let the finger vent before the run
            plant.advance(1 / cfg.loop.rate_hz)
        finger.estimator.update(plant.read_pressure())

    t0 = plant.now()
    schedule = [(t0 + i * T_STEP_S, k) for i, k in enumerate(a.targets)]
    t_end = t0 + len(a.targets) * T_STEP_S
    log = run_closed_loop(plant, finger, schedule, t_end, cfg.loop.rate_hz)
    t = log["t"]

    print(f"\n{'step':>4} {'kappa*':>7} {'P*':>6} | {'settle s':>7} {'over kPa':>7} {'err kPa':>7}"
          f" | {'settle s':>7} {'over 1/m':>7} {'err 1/m':>7}")
    prev_p = 0.0
    for i, (ts, kap) in enumerate(schedule):
        t_next = t_end if i + 1 == len(schedule) else schedule[i + 1][0]
        Pt = finger.map.pressure_for(finger.mapping.target_curvature(kap))
        Kt = float(finger.map.kappa(Pt))
        if Pt == prev_p:
            continue
        mp = step_metrics(t, log["P_true"], ts, t_next, Pt, P_BAND_KPA)
        mk = step_metrics(t, log["kappa_true"], ts, t_next, Kt, K_BAND)
        print(f"{i:>4} {kap:7.1f} {Pt:6.2f} | {fmt(mp)} | {fmt(mk)}")
        prev_p = Pt
    print(f"\npeak true pressure {log['P_true'].max():.2f} kPa (cap {cfg.safety.p_cap_kPa}); "
          f"settle band +-{P_BAND_KPA} kPa / +-{K_BAND} 1/m, on true (not measured) values")

    if not a.no_plot:
        plot(log, t0, a.out, cfg.controller.type)


def fmt(m):
    settle = "never" if m["settle_s"] is None else f"{m['settle_s']:.2f}"
    return f"{settle:>7} {m['overshoot']:7.2f} {m['steady_err']:+7.2f}"


def plot(log, t0, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = log["t"] - t0
    fig, ax = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    ax[0].plot(t, log["P_true"], label="P (true)")
    ax[0].plot(t, log["P_meas"], ".", ms=2, alpha=0.5, label="P_meas (calibrated)")
    ax[0].plot(t, log["P_target"], "k--", label="P*")
    ax[0].set_ylabel("pressure (kPa)")
    ax[1].plot(t, log["kappa_true"], label="kappa (true)")
    ax[1].plot(t, log["kappa_target"], "k--", label="kappa*")
    ax[1].set_ylabel("curvature (1/m)")
    ax[2].plot(t, log["duty"], label="duty")
    ax[2].set_ylabel("duty")
    ax[2].set_xlabel("time (s)")
    for x in ax:
        x.legend(loc="upper right", fontsize=8)
        x.grid(alpha=0.3)
    ax[0].set_title(f"Finger pressure control: {title}")
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"saved {path}")


if __name__ == "__main__":
    main()
