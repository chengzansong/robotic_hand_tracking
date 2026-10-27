# Control

Turns the curvature perception reports for each hand finger into a pump duty for the matching robot finger. First working version: it drives a simulated finger (`LumpedPlant`) directly, on the PC. There is no serial link or firmware protocol yet.

```
perception curvature (1/m)
  -> mapping.py             human curvature -> robot target curvature kappa*
  -> shape/pressure_map.py  kappa* -> pressure P*, from the FEM table, clipped to the pressure cap
  -> controllers/           (P*, P_hat, dP_hat) -> duty u      pi_ff (default), bang_bang, open_loop_volume
  -> controllers/safety.py  pressure cap, max duty, watchdog, u = 0 on a fault
  -> Plant.set_duty()
Plant.read_pressure() -> estimators/  calibration, median-of-3, alpha-beta filter, fault flags -> P_hat, dP_hat
```

The design and the operator's decisions are in [plans/2026-10-26_control-plan.md](../plans/2026-10-26_control-plan.md) (section 8). The data formats are in [CONTRACTS.md](../CONTRACTS.md).

## Layout

| Path | Contents |
| --- | --- |
| `core/` | `plant.py`: the `Plant` protocol. `config.py`: loads the TOML. `types.py`: `Sample` and `Fault`. `pneumatics.py`: pump and bleed flow laws |
| `plants/lumped.py` | `LumpedPlant`: pump, finger and fixed bleed as one ODE, plus the sensor model |
| `estimators/` | `calibration.py` (gain, offset, startup zero), `filters.py`, `pressure.py` (the estimator and fault flags), `leak_id.py` (fits the bleed law from a decay) |
| `controllers/` | `pi_ff.py` (default; `pi_ff_valve` is a stub), `bang_bang.py`, `open_loop_volume.py`, `safety.py` |
| `shape/` | `export_map.py` (reads a Stage 1 `results.json`), `pressure_map.py` (interpolation and inverse), `maps/` (committed tables) |
| `mapping.py` | Perception curvature to robot target |
| `finger.py` | `FingerController`: one finger's whole stack behind `tick(t, samples, kappa_h) -> duty` |
| `runner.py` | Closed-loop runner and step-response metrics |
| `configs/finger_eco30.toml` | Every tunable in one commented file |
| `demo.py`, `tests/` | The demo and the tests |

Each finger picks its controller with `type` under `[controller]` in its TOML.

## Run

From the repo root, with the project's `.venv`:

```
.venv\Scripts\python -m control.demo
.venv\Scripts\python -m control.demo --controller bang_bang --bleed-mm 0.3
.venv\Scripts\python -m control.demo --identify-leak --bleed-mm 0.4
.venv\Scripts\python -m pytest control/tests -q
```

The demo zeroes the sensor, then steps one finger through curvature targets of 0, 17, 38, 27 and 0 1/m (`--targets` changes them). It prints settle time, overshoot and steady error for each step, in pressure and in curvature, judged on the plant's true values and not the noisy sensor. It saves a plot (pressure, measured pressure, target, curvature, duty) to `control/outputs/demo.png`, which git ignores. `--bleed-mm` changes the plant's bleed only. The controller keeps its belief from `[model]`, so a different size shows how it copes with a wrong leak model. `--identify-leak` first runs a decay test and gives the controller the fitted leak.

On the default settings, `pi_ff` settles in about 0.7 to 1.4 s per step with under 0.2 kPa overshoot. Releasing to 0 takes about 3 s, because the bleed is the only way out. `bang_bang` ripples by about 0.3 kPa around the target, so it rarely sits inside a tight band.

## Shape map

`shape/maps/finger_poc_eco30_meshB.json` is the pressure-to-shape table for the `finger_poc` design in Ecoflex 00-30. It comes from the Stage 1 free sweep `simulation/outputs/eco30_sweep_10kPa/20261020T121529Z/results.json`. Control reads only this JSON file and never imports `simulation/`. To make a table for another design or material:

```
.venv\Scripts\python -m control.shape.export_map <results.json> -o control/shape/maps/<name>.json --design <id> --material <id>
```

The table holds pressure, bend angle, curvature, cavity volume and blocked force (null if the run has none). Curvature is 1/R_fit from the run's `turn_radius` rows once the span turns past about 90 degrees, 1/R_avg below that, and bend angle divided by arc length for any level with no row. The `curvature_source` list in the file records the rule used at each point. The source run's `turn_radius` rows give 1/R_avg at 0.5 and 1 kPa and 1/R_fit from 2 to 5 kPa (17.1, 27.0, 38.4 and 54.5 1/m). For older runs that have no `turn_radius` rows, `--radius-mm 2:58.4,3:37.1` supplies the radii by hand. The cavity volume gives the finger's compliance for the plant. The source is a free sweep, so the table has no blocked force.

The map never extrapolates. The target pressure is clipped to the lower of the table's top level (5 kPa, 54 1/m) and the design cap in the config.

## Placeholders

Nothing here is measured on the real pump or sensor yet. Each constant is in `configs/finger_eco30.toml` and is replaced by a bench test.

| Constant | Value | Replaced by |
| --- | --- | --- |
| Pump free flow Q0 | 8 ml/s (the only bench point is 5.0 ml/s at an unknown back-pressure) | Bench B2-B4 |
| Pump stall pressure | 30 kPa | B2-B4 |
| Pump start duty u_min | 0.25 | B2-B4 |
| Pump spin-up lag | 0.05 s (est.) | B2-B4 |
| Bleed size and law | 0.3 mm orifice, n = 0.5, Cd = 0.6 | B5, or the leak identification routine |
| Tube volume | 1 ml | Geometry as built |
| Sensor gain, offset, noise | 1% gain error, 0.15 kPa offset, 0.08 kPa noise, 0.01 kPa LSB, 30 ms latency | B6 |
| Estimator gain | 1.0 | B6 calibration against a reference gauge |
| Fault thresholds | stuck after 20 identical readings, 25 kPa/s slope, range -1 to 12 kPa | B6 (the stuck rule assumes noise above one LSB) |
| Controller model of pump and leak | Same numbers as the plant | B2-B5; the loop also learns u_ff |
| Controller gains | Kp 0.7, Ki 0.5, Kd 0.05 | Tuned on `LumpedPlant`; retune in the real-time sim (phase C4) and on the bench |
| Pressure cap | 5.5 kPa | FEM and the bench, per material |
| Human finger length | 85 mm, used in angle mode only | Measured per user |

The finger's shape and compliance come from the FEM, so the bench B7 will correct them. The shape map is quasi-static: the finger has no bending lag in `LumpedPlant`.

With the placeholder pump, a 0.4 mm bleed cannot hold pressures above about 3.5 kPa at the 0.9 duty cap. The pump runs out of flow, so the loop saturates and the demo reports "never" for the 38 1/m step. A larger pump or a smaller bleed fixes it.

## Plug in another plant

Control talks to a plant only through the protocol in `core/plant.py`:

```python
class Plant(Protocol):
    def set_duty(self, finger: str, duty: float) -> None: ...   # 0-1, held until changed
    def read_pressure(self) -> list[Sample]: ...                # new raw samples, oldest first
    def now(self) -> float: ...                                 # plant time, s
    def advance(self, dt_s: float) -> None: ...                 # lockstep plants only
```

A `Sample` is `(t_s, finger, kPa, seq)`: the raw reading, the time it was taken, and a sequence number. The real-time simulator in `simulation/stage2_hand/` (`SimIO`) and the future ESP32 bridge should implement the same calls. Control never imports `simulation/`; the simulator may import `control.core`. To drive another plant, build a `FingerController` and call it each control tick:

```python
finger = FingerController.from_toml("control/configs/finger_eco30.toml")
duty = finger.tick(plant.now(), plant.read_pressure(), kappa_h)    # kappa_h: perception curvature, 1/m
plant.set_duty(finger.name, duty)
```

`runner.run_closed_loop` does this at a fixed rate and logs the result. A real-time plant makes `advance` a no-op or a sleep. The ground-truth columns of the log are filled only if the plant has `true_pressure` and `true_curvature`.

## Not built yet

- `SerialESP32`: no protocol is fixed. The PC-to-microcontroller link is a later step.
- `pi_ff_valve`: a stub, because v1 has no release valve.
- The vision trim, which would adjust P* until the observed robot curvature matches the target.
- A maximum continuous on-time for the pump, and a multi-finger loop (each finger has its own `FingerController` and TOML).
