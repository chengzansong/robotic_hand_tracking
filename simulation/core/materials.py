"""Load material cards from simulation/materials/.

Cards are solver independent. Backends convert the canonical parameters
documented in simulation/materials/README.md into their own conventions.

Two kinds of card:
- kind = "solid" (the default): a bulk hyperelastic material (silicones).
- kind = "sheet": a thin reinforcement such as a woven fabric, described by
  its in-plane tension per unit width. It is not a solid on its own; a study
  casts it into a part made of a matrix solid (the silicone that impregnates
  it). Model "fabric_pow_linear", canonical form (per yarn family, warp and
  weft, stretch lam along the yarn, T = nominal tension per unit width):
      T = 0                                for lam <= 1 (buckles: no compression)
      toe, 1 < lam < lam0:   T = S c_t lam (lam^2 - 1)^(beta - 1)
      linear, lam >= lam0:   dT/dlam = S, continuous in value and slope at lam0
  with S the linear-region stiffness per width (N/mm), lam0 the toe stretch
  (yarn crimp), and beta >= 2 the toe exponent; c_t follows from slope
  continuity. This is FEBio's "fiber-pow-linear" scaled by the modelled layer
  thickness (backends/febio/materials.py derives c_t and the FEBio modulus).
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
# Sheet (fabric) models: required keys; optional keys are listed in SHEET_OPTIONAL
SHEET_MODEL_KEYS = {
    "fabric_pow_linear": ("stiffness_N_per_mm", "toe_stretch", "toe_exponent"),
}
SHEET_OPTIONAL = ("stiffness_weft_N_per_mm", "thickness_mm", "areal_mass_g_m2")
KINDS = ("solid", "sheet")


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

    @property
    def is_sheet(self):
        return self.model in SHEET_MODEL_KEYS

    def warp_weft_N_per_mm(self):
        """Sheet fits: (warp, weft) linear-region tension stiffness per width."""
        if not self.is_sheet:
            raise ValueError(f"{self.model} is not a sheet model")
        s = self.params["stiffness_N_per_mm"]
        return float(s), float(self.params.get("stiffness_weft_N_per_mm", s))

    def initial_shear_modulus_kPa(self):
        """Small-strain shear modulus implied by the fit."""
        if self.is_sheet:
            raise ValueError(f"{self.model} is a sheet model; it has no shear modulus")
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
    kind: str = "solid"

    @property
    def default_fit(self):
        return self.fits[0]


def _parse_fit(raw, card_id, kind="solid"):
    model = raw.get("model")
    keys = SHEET_MODEL_KEYS if kind == "sheet" else MODEL_KEYS
    if model not in keys:
        raise ValueError(f"{card_id}: unknown {kind} model {model!r}")

    missing = [k for k in keys[model] if k not in raw]
    if missing:
        raise ValueError(f"{card_id}: {model} fit is missing {missing}")

    params = {k: raw[k] for k in keys[model]}
    if kind == "sheet":
        params.update({k: raw[k] for k in SHEET_OPTIONAL if k in raw})
        if any(float(raw[k]) <= 0 for k in ("stiffness_N_per_mm", "stiffness_weft_N_per_mm")
               if k in raw):
            raise ValueError(f"{card_id}: sheet stiffness must be positive")
        if not float(params["toe_stretch"]) > 1.0:
            raise ValueError(f"{card_id}: toe_stretch must be > 1")
        if float(params["toe_exponent"]) < 2.0:
            raise ValueError(f"{card_id}: toe_exponent must be >= 2 (FEBio fiber-pow-linear)")
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
    kind = raw.get("kind", "solid")
    if kind not in KINDS:
        raise ValueError(f"{path.name}: kind must be one of {KINDS}, got {kind!r}")

    return Material(
        id=raw["id"],
        name=raw.get("name", raw["id"]),
        family=raw.get("family", ""),
        density_kg_m3=raw.get("density_kg_m3", 1000.0),
        fits=[_parse_fit(f, material_id, kind) for f in raw["fits"]],
        datasheet=raw.get("datasheet", {}),
        kind=kind,
    )


def card_kind(material_id, materials_dir=MATERIALS_DIR):
    with open(Path(materials_dir) / f"{material_id}.toml", "rb") as f:
        return tomllib.load(f).get("kind", "solid")


def list_materials(materials_dir=MATERIALS_DIR, include_test=False, kind="solid"):
    """Ids of the cards in the materials directory of one kind ("solid" by
    default, "sheet", or None for every card)."""
    ids = sorted(p.stem for p in Path(materials_dir).glob("*.toml"))
    ids = ids if include_test else [i for i in ids if not i.startswith("_")]
    return ids if kind is None else [i for i in ids if card_kind(i, materials_dir) == kind]
