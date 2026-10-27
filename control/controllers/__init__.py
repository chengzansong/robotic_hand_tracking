"""L2 pressure controllers and the safety layer.

Every controller has `step(t, P_hat, dP_hat, P_target) -> u` and `reset()`.
Pick one per finger with `type` in the TOML, and build it with `build_controller`.
"""

from .bang_bang import BangBang
from .open_loop_volume import OpenLoopVolume
from .pi_ff import PIFF, PIFFValve, model_ff_table
from .safety import Safety

CONTROLLERS = {
    "bang_bang": BangBang,
    "open_loop_volume": OpenLoopVolume,
    "pi_ff": PIFF,
    "pi_ff_valve": PIFFValve,
}


def build_controller(cfg, shape_map):
    """Build the controller named in `cfg.controller` from its [controller.<type>] table."""
    kind = cfg.controller.type
    if kind not in CONTROLLERS:
        raise ValueError(f"unknown controller {kind!r}; choose from {sorted(CONTROLLERS)}")
    return CONTROLLERS[kind](shape_map=shape_map, model=cfg.model,
                             p_max_kPa=shape_map.p_max, **cfg.controller.params)


__all__ = ["BangBang", "OpenLoopVolume", "PIFF", "PIFFValve", "Safety",
           "CONTROLLERS", "build_controller", "model_ff_table"]
