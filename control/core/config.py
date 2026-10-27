"""Load a finger's TOML config into plain dataclasses.

One file holds every tunable (see control/configs/finger_eco30.toml). Unknown
keys are an error, so a typo in the TOML cannot silently fall back to a default.
"""

import tomllib
from dataclasses import dataclass, field, fields
from pathlib import Path

CONTROL_DIR = Path(__file__).resolve().parent.parent


def _make(cls, d, where):
    names = {f.name for f in fields(cls)}
    bad = set(d) - names
    if bad:
        raise ValueError(f"unknown key(s) in [{where}]: {sorted(bad)}")
    return cls(**d)


@dataclass
class FingerCfg:
    name: str = "Index"
    shape_map: str = "shape/maps/finger_poc_eco30_meshB.json"
    p_cap_kPa: float = 5.5


@dataclass
class MappingCfg:
    mode: str = "curvature"       # "curvature" (direct) or "angle"
    L_h_mm: float = 85.0
    L_robot_mm: float = 112.0


@dataclass
class LoopCfg:
    rate_hz: float = 40.0


@dataclass
class PumpCfg:
    q0_ml_s: float = 8.0
    p_stall_kPa: float = 30.0
    u_min: float = 0.25
    tau_motor_s: float = 0.05


@dataclass
class BleedCfg:
    diameter_mm: float = 0.3
    cd: float = 0.6
    n: float = 0.5
    k_ml_s: float | None = None   # overrides diameter_mm when set (ml/s per kPa^n)


@dataclass
class PlantCfg:
    v_tube_ml: float = 1.0
    dt_s: float = 0.001
    watchdog_s: float = 0.2
    seed: int = 1


@dataclass
class SensorCfg:
    rate_hz: float = 40.0
    latency_s: float = 0.03
    gain_error_sigma: float = 0.01
    offset_sigma_kPa: float = 0.15
    noise_sigma_kPa: float = 0.08
    lsb_kPa: float = 0.01
    range_kPa: float = 40.0


@dataclass
class EstimatorCfg:
    zero_s: float = 1.0
    gain: float = 1.0
    offset_kPa: float = 0.0
    filter: str = "alpha_beta"    # "alpha_beta" or "ema"
    alpha: float = 0.3
    beta: float = 0.05
    ema_alpha: float = 0.4
    range_lo_kPa: float = -1.0
    range_hi_kPa: float = 12.0
    stuck_n: int = 20
    max_slope_kPa_s: float = 25.0


@dataclass
class ModelCfg:
    """What the controller believes about the pump and leak (its prior)."""

    q0_ml_s: float = 8.0
    p_stall_kPa: float = 30.0
    u_min: float = 0.25
    bleed_mm: float = 0.3
    cd: float = 0.6
    n: float = 0.5
    k_ml_s: float | None = None   # overrides bleed_mm when set (e.g. from leak identification)


@dataclass
class SafetyCfg:
    p_cap_kPa: float = 5.5
    cap_release_kPa: float = 0.3
    max_duty: float = 0.9
    watchdog_s: float = 0.2


@dataclass
class ControllerCfg:
    """`type` picks the algorithm; `tables` holds every [controller.<type>] table."""

    type: str = "pi_ff"
    tables: dict = field(default_factory=dict)

    @property
    def params(self):
        return self.tables.get(self.type, {})


@dataclass
class Config:
    finger: FingerCfg
    mapping: MappingCfg
    loop: LoopCfg
    pump: PumpCfg
    bleed: BleedCfg
    plant: PlantCfg
    sensor: SensorCfg
    estimator: EstimatorCfg
    model: ModelCfg
    safety: SafetyCfg
    controller: ControllerCfg
    path: Path | None = None

    def shape_map_path(self):
        p = Path(self.finger.shape_map)
        return p if p.is_absolute() else CONTROL_DIR / p


_SECTIONS = {
    "finger": FingerCfg, "mapping": MappingCfg, "loop": LoopCfg, "pump": PumpCfg,
    "bleed": BleedCfg, "plant": PlantCfg, "sensor": SensorCfg,
    "estimator": EstimatorCfg, "model": ModelCfg, "safety": SafetyCfg,
}


def config_from_dict(raw, path=None):
    raw = dict(raw)
    kw = {name: _make(cls, raw.pop(name, {}), name) for name, cls in _SECTIONS.items()}
    ctl = dict(raw.pop("controller", {}))
    kw["controller"] = ControllerCfg(type=ctl.pop("type", "pi_ff"), tables=ctl)
    if raw:
        raise ValueError(f"unknown section(s): {sorted(raw)}")
    return Config(path=path, **kw)


def load_config(path):
    """Read a finger TOML file into a Config."""
    path = Path(path)
    with open(path, "rb") as f:
        return config_from_dict(tomllib.load(f), path=path)
