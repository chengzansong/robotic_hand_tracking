"""Load material cards from simulation/materials/.

Cards are solver independent. Backends convert the canonical parameters
documented in simulation/materials/README.md into their own conventions.
"""

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

MATERIALS_DIR = Path(__file__).resolve().parent.parent / "materials"

# Required parameter keys for each constitutive model, all in kPa
MODEL_KEYS = {
    "neo_hookean": ("mu_kPa",),
    "mooney_rivlin": ("c10_kPa", "c01_kPa"),
    "yeoh": ("c10_kPa", "c20_kPa", "c30_kPa"),
    "ogden": ("mu_kPa", "alpha"),
}


@dataclass
class Fit:
    model: str
    params: dict
    source: str
    confidence: str
    bulk_modulus_kPa: float | None = None
    poisson: float | None = None
    stretch_range: tuple[float, float] | None = None
    notes: str = ""

    def initial_shear_modulus_kPa(self):
        """Small-strain shear modulus implied by the fit."""
        p = self.params
        if self.model == "neo_hookean":
            return p["mu_kPa"]
        if self.model == "mooney_rivlin":
            return 2 * (p["c10_kPa"] + p["c01_kPa"])
        if self.model == "yeoh":
            return 2 * p["c10_kPa"]
        if self.model == "ogden":
            return sum(p["mu_kPa"])
        raise ValueError(f"unknown model {self.model}")


@dataclass
class Material:
    id: str
    name: str
    family: str
    density_kg_m3: float
    fits: list[Fit]
    datasheet: dict = field(default_factory=dict)

    @property
    def default_fit(self):
        return self.fits[0]


def _parse_fit(raw, card_id):
    model = raw.get("model")
    if model not in MODEL_KEYS:
        raise ValueError(f"{card_id}: unknown model {model!r}")

    missing = [k for k in MODEL_KEYS[model] if k not in raw]
    if missing:
        raise ValueError(f"{card_id}: {model} fit is missing {missing}")

    params = {k: raw[k] for k in MODEL_KEYS[model]}
    if model == "ogden" and len(params["mu_kPa"]) != len(params["alpha"]):
        raise ValueError(f"{card_id}: ogden mu_kPa and alpha differ in length")

    stretch = raw.get("stretch_range")
    return Fit(
        model=model,
        params=params,
        source=raw.get("source", ""),
        confidence=raw.get("confidence", "low"),
        bulk_modulus_kPa=raw.get("bulk_modulus_kPa"),
        poisson=raw.get("poisson"),
        stretch_range=tuple(stretch) if stretch else None,
        notes=raw.get("notes", ""),
    )


def load_material(material_id, materials_dir=MATERIALS_DIR):
    """Load one card by id, for example "ecoflex_00-30"."""
    path = Path(materials_dir) / f"{material_id}.toml"
    with open(path, "rb") as f:
        raw = tomllib.load(f)

    if raw.get("id") != material_id:
        raise ValueError(f"{path.name}: id {raw.get('id')!r} does not match the file name")
    if not raw.get("fits"):
        raise ValueError(f"{path.name}: no [[fits]] entries")

    return Material(
        id=raw["id"],
        name=raw.get("name", raw["id"]),
        family=raw.get("family", ""),
        density_kg_m3=raw.get("density_kg_m3", 1000.0),
        fits=[_parse_fit(f, material_id) for f in raw["fits"]],
        datasheet=raw.get("datasheet", {}),
    )


def list_materials(materials_dir=MATERIALS_DIR, include_test=False):
    """Ids of every card in the materials directory."""
    ids = sorted(p.stem for p in Path(materials_dir).glob("*.toml"))
    return ids if include_test else [i for i in ids if not i.startswith("_")]
