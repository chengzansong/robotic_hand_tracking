"""compare_gif: level matching, blocked-force lookup, and frame annotation (no rendering)."""

import numpy as np
import pytest

from simulation.stage1_finger.compare_gif import annotate, blocked_force_at, common_levels


def test_common_levels_and_cap():
    a = {"p": np.array([0, 0.5, 1, 2, 3, 4, 5, 6])}
    b = {"p": np.array([0, 1, 2, 3, 4, 5])}
    assert common_levels(a, b) == [0, 1, 2, 3, 4, 5]
    assert common_levels(a, b, max_kPa=3) == [0, 1, 2, 3]


def test_blocked_force_lookup():
    blk = {"p": np.array([0.0, 1, 2, 3, 4, 5]), "F": np.array([0, 2.8, 6.3, 9.9, 13.6, 17.3]) / 1000,
           "mesh_name": "mesh A"}
    f, exact = blocked_force_at(blk, 3.0)
    assert f == pytest.approx(0.0099) and exact
    f, exact = blocked_force_at(blk, 0.5)
    assert f == pytest.approx(0.0014) and not exact
    assert blocked_force_at(blk, 6.0) == (None, False)
    assert blocked_force_at(None, 1.0) == (None, False)


def test_annotation_text():
    s = {"pressure_kPa": 4.0, "bend_deg": 208.6, "R_mm": 26.0, "R_kind": "fit",
         "kappa_per_m": 1000 / 26.0, "underside_length": 0.934}
    blk = {"p": np.array([0.0, 4.0]), "F": np.array([0.0, 0.0136]), "mesh_name": "mesh A"}
    t = annotate("Baseline", s, 0.0136, True, blk, "mesh B")
    assert "p = 4 kPa" in t and "bend 209 deg" in t
    assert "26.0 mm (fit)" in t and "38.5 1/m" in t
    assert "13.6 mN, mesh A" in t
    t0 = annotate("B", {**s, "pressure_kPa": 0.0, "bend_deg": 0.0, "R_mm": None}, 0.0, True, blk,
                  "mesh A")
    assert "curvature 0 1/m" in t0 and "0.0 mN" in t0 and "mesh A" not in t0.split("mN")[1]
    assert "not run" in annotate("B", s, None, False, None, "mesh B")
