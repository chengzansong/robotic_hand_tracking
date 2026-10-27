"""Shape map, its inverse, the exporter and the perception mapping."""

import json
import math

import numpy as np
import pytest

from control.finger import FingerController
from control.mapping import Mapping
from control.shape import export_map
from control.shape.pressure_map import ShapeMap


@pytest.fixture
def smap(cfg):
    return ShapeMap.load(cfg.shape_map_path(), cfg.finger.p_cap_kPa)


def test_table_matches_the_fem_points(smap):
    for p, k in [(2, 17.1), (3, 27.0), (4, 38.5), (5, 54.3)]:
        assert smap.kappa(p) == pytest.approx(k, abs=0.2)
    assert smap.bend_deg(4.0) == pytest.approx(208.6, abs=0.1)


def test_inverse_round_trips(smap):
    for P in np.linspace(0.0, smap.p_max, 40):
        assert smap.pressure_for(smap.kappa(P)) == pytest.approx(P, abs=1e-3)
    for k in np.linspace(0.0, smap.kappa_max, 40):
        assert float(smap.kappa(smap.pressure_for(k))) == pytest.approx(k, abs=1e-2)


def test_inverse_is_monotone_and_clipped(smap, cfg):
    ks = np.linspace(-5, 200, 300)
    ps = [smap.pressure_for(k) for k in ks]
    assert np.all(np.diff(ps) >= -1e-12)
    assert ps[0] == 0.0
    assert max(ps) <= min(cfg.finger.p_cap_kPa, smap.p[-1]) + 1e-12


def test_cap_below_table_top_clips_pressure(cfg):
    capped = ShapeMap.load(cfg.shape_map_path(), p_cap_kPa=3.0)
    assert capped.pressure_for(100.0) == pytest.approx(3.0)


def test_compliance_is_in_the_expected_range(smap):
    c = [smap.compliance_ml_per_kPa(P) for P in (0.5, 2.0, 4.0)]
    assert all(1.0 < x < 6.0 for x in c)         # control plan: 1.7-4.3 ml/kPa


def test_exporter_prefers_turn_radius_and_falls_back(tmp_path):
    steps = [{"pressure_kPa": p, "bend_angle_deg": b, "cavity_volume_mm3": 1000.0 + 100 * p,
              "blocked_force_N": None, "requested_level": True}
             for p, b in [(0.0, 0.0), (1.0, 50.0), (2.0, 100.0), (3.0, 150.0)]]
    turn = [{"pressure_kPa": 2.0, "R_fit_mm": 50.0, "R_avg_mm": 60.0, "fit_reliable": True},
            {"pressure_kPa": 1.0, "R_fit_mm": 90.0, "R_avg_mm": 100.0, "fit_reliable": False}]
    path = tmp_path / "results.json"
    path.write_text(json.dumps({"design": "d", "steps": steps, "turn_radius": turn}))
    t = export_map.export_map(path, material="m")
    assert t["curvature_per_m"][:3] == [0.0, 10.0, 20.0]       # 1/R_avg, then 1/R_fit
    assert t["curvature_source"] == ["zero", "1/R_avg", "1/R_fit",
                                     f"bend/arc_length ({t['arc_length_mm']:.1f} mm)"]
    expect = math.radians(150) / (t["arc_length_mm"] / 1000)
    assert t["curvature_per_m"][3] == pytest.approx(expect, abs=1e-3)
    assert t["cavity_volume_mm3"][0] == 1000.0 and t["blocked_force_N"] is None


def test_mapping_direct_is_identity():
    m = Mapping()
    assert m.target_curvature(17.0) == 17.0
    assert m.target_curvature(-3.0) == 0.0


def test_mapping_angle_mode_keeps_the_angle():
    m = Mapping("angle", L_h_mm=85.0, L_robot_mm=112.0)
    k = m.target_curvature(30.0)
    assert k == pytest.approx(30.0 * 85.0 / 112.0)
    assert k * 0.112 == pytest.approx(30.0 * 0.085)       # same bend angle in radians
    assert k < 30.0                                       # robot is longer, so lower curvature


def test_mapping_rejects_unknown_mode():
    with pytest.raises(ValueError):
        Mapping("nonsense")


def test_angle_mode_gives_less_pressure_end_to_end(cfg):
    direct = FingerController(cfg)
    cfg.mapping.mode = "angle"
    angle = FingerController(cfg)
    direct.tick(0.0, [], 30.0)
    angle.tick(0.0, [], 30.0)
    assert angle.P_target < direct.P_target
