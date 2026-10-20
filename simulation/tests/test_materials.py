"""Material card -> FEBio conversion keeps the initial shear modulus."""

import pytest

from simulation.core.materials import Fit, list_materials, load_material
from simulation.stage1_finger.backends.febio.materials import (
    UnsupportedMaterial, bulk_modulus_kPa, febio_initial_shear_MPa, febio_material)

SYNTHETIC = [
    Fit("neo_hookean", {"mu_kPa": 40.0}, "", "low"),
    Fit("mooney_rivlin", {"c10_kPa": 15.0, "c01_kPa": 5.0}, "", "low"),
    Fit("yeoh", {"c10_kPa": 12.58, "c20_kPa": 0.5759, "c30_kPa": 0.01}, "", "low"),
    Fit("ogden", {"mu_kPa": [14.55], "alpha": [2.34]}, "", "low"),
    Fit("ogden", {"mu_kPa": [20.0, -1.0, 0.5], "alpha": [1.3, -2.0, 5.0]}, "", "low"),
]


def _check(fit):
    mtype, params, _ = febio_material(fit, bulk_to_shear=100.0)
    mu0_MPa = febio_initial_shear_MPa(mtype, params)
    assert mu0_MPa == pytest.approx(fit.initial_shear_modulus_kPa() * 1e-3, rel=1e-12)
    assert params["k"] > 0
    return mtype, params


@pytest.mark.parametrize("fit", SYNTHETIC, ids=lambda f: f.model)
def test_initial_shear_modulus_preserved(fit):
    _check(fit)


def test_every_card_and_fit_converts():
    ids = list_materials(include_test=True)
    assert "_placeholder_neo_hookean" in ids
    for material_id in ids:
        for fit in load_material(material_id).fits:
            _check(fit)


def test_febio_keywords():
    assert _check(SYNTHETIC[0]) == ("Mooney-Rivlin", pytest.approx(
        {"c1": 0.02, "c2": 0.0, "k": 4.0}))
    mtype, params = _check(SYNTHETIC[2])
    assert mtype == "Yeoh" and params["c1"] == pytest.approx(0.01258)
    mtype, params = _check(SYNTHETIC[3])
    assert mtype == "Ogden"
    assert params["c1"] == pytest.approx(2 * 0.01455) and params["m1"] == 2.34


def test_bulk_modulus_sources():
    nh = Fit("neo_hookean", {"mu_kPa": 40.0}, "", "low")
    assert bulk_modulus_kPa(nh, 50.0) == (2000.0, "study bulk_to_shear 50.0")
    nh.poisson = 0.49
    k, src = bulk_modulus_kPa(nh, 50.0)
    assert k == pytest.approx(2 * 40 * 1.49 / (3 * 0.02)) and "poisson" in src
    nh.bulk_modulus_kPa = 1234.0
    assert bulk_modulus_kPa(nh, 50.0)[0] == 1234.0


def test_unsupported_raises():
    with pytest.raises(UnsupportedMaterial):
        febio_material(Fit("ogden", {"mu_kPa": [1.0] * 7, "alpha": [2.0] * 7}, "", "low"), 100)
    with pytest.raises(UnsupportedMaterial):
        febio_material(Fit("arruda_boyce", {}, "", "low"), 100)
