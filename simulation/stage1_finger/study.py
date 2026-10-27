"""Load a Stage 1 study TOML: design + material fit per part + mesh + load case.

Format (see studies/*.toml):

    name = "poc_smoke"                 # optional, defaults to the file name
    design = "finger_poc"              # folder in simulation/models/
    backend = "febio"                  # key of backends.BACKENDS
    bulk_to_shear = 100.0              # K/mu0, used only when a fit gives neither
                                       # bulk_modulus_kPa nor poisson

    [materials]                        # per-part override of the design defaults
    top = "ecoflex_00-30"              # fit 0 (the card's default)
    base = { material = "ecoflex_00-30", fit = 2 }   # fit by index
    # A sheet card (kind = "sheet", e.g. a fabric) on a layer part such as
    # "cloth" is cast into a matrix silicone: by default the base part's
    # material; stiffness_scale multiplies the sheet's tension stiffness.
    cloth = { material = "fabric_woven_cotton", matrix = "ecoflex_00-30", stiffness_scale = 1.0 }

    [mesh]
    size_mm = 2.0
    wall_size_mm = 0.8
    wall_distance_mm = 2.0

    symmetry = ""                      # "" or "half_y" (half model, y mid-plane)

    [load]
    case = "free"                      # free | blocked_force
    pressure_kPa = 20.0                # max chamber gauge pressure
    steps = 10                         # output points, evenly spaced in pressure
    # or explicit output levels (then pressure_kPa and steps are ignored):
    # levels_kPa = [0.5, 1, 2, 3, 5, 7.5, 10, 15, 20]

    [contact]                          # slot-wall contact; on by default
    enabled = true
    fric_coeff = 0.0                   # see backends/febio/writer.py CONTACT_DEFAULTS

    [solver]                           # optional backend settings, see the backend
"""

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

from simulation.core.designs import load_design
from simulation.core.materials import MATERIALS_DIR, list_materials, load_material

from .mesh import MeshSpec

LOAD_CASES = ("free", "blocked_force")


class StudyError(ValueError):
    pass


@dataclass
class PartMaterial:
    material: object       # simulation.core.materials.Material
    fit_index: int
    # sheet materials only: the solid that impregnates the sheet, a scale on the
    # sheet's tension stiffness, and the layer geometry from the design part
    # (thickness_mm, warp_axis, weft_axis)
    matrix: "PartMaterial | None" = None
    stiffness_scale: float = 1.0
    layer: dict | None = None

    @property
    def fit(self):
        return self.material.fits[self.fit_index]

    @property
    def is_sheet(self):
        return self.material.kind == "sheet"

    def describe(self):
        f = self.fit
        out = {"material": self.material.id, "fit_index": self.fit_index,
               "model": f.model, "params": f.params, "source": f.source,
               "confidence": f.confidence, "bulk_modulus_kPa": f.bulk_modulus_kPa,
               "poisson": f.poisson, "stretch_range": f.stretch_range,
               "initial_shear_modulus_kPa": None if self.is_sheet else f.initial_shear_modulus_kPa(),
               "density_kg_m3": self.material.density_kg_m3}
        if self.is_sheet:
            out.update(kind="sheet", stiffness_scale=self.stiffness_scale, layer=self.layer,
                       matrix=self.matrix.describe())
        return out


@dataclass
class LoadCase:
    case: str = "free"
    pressure_kPa: float = 10.0
    steps: int = 10
    levels_kPa: list | None = None

    def levels(self):
        """Output pressure levels in kPa, ascending."""
        if self.levels_kPa:
            return [float(p) for p in self.levels_kPa]
        return [self.pressure_kPa * k / self.steps for k in range(1, self.steps + 1)]

    @property
    def max_kPa(self):
        return self.levels()[-1]

    def fractions(self):
        """Load fractions t in (0, 1] of each level; pressure is max_kPa * t."""
        return [p / self.max_kPa for p in self.levels()]


@dataclass
class Study:
    name: str
    path: Path
    design: object         # simulation.core.designs.Design
    materials: dict        # part name -> PartMaterial
    mesh: MeshSpec
    load: LoadCase
    backend: str = "febio"
    bulk_to_shear: float = 100.0
    solver: dict = field(default_factory=dict)
    contact: dict = field(default_factory=dict)
    raw: dict = field(default_factory=dict)


SHEET_KEYS = ("material", "fit", "matrix", "stiffness_scale")


def _material(part, entry, materials_dir):
    if isinstance(entry, str):
        material_id, fit_index = entry, 0
    elif isinstance(entry, dict) and "material" in entry:
        material_id, fit_index = entry["material"], int(entry.get("fit", 0))
        unknown = set(entry) - set(SHEET_KEYS)
        if unknown:
            raise StudyError(f"[materials] {part}: unknown keys {sorted(unknown)}")
    else:
        raise StudyError(f"[materials] {part}: use a card id or "
                         f"{{ material = \"id\", fit = 0 }}, got {entry!r}")
    if not (Path(materials_dir) / f"{material_id}.toml").exists():
        raise StudyError(
            f"material card {material_id!r} for part {part!r} not found in {materials_dir}. "
            f"Available: {', '.join(list_materials(materials_dir, include_test=True, kind=None))}")
    material = load_material(material_id, materials_dir)
    if not 0 <= fit_index < len(material.fits):
        raise StudyError(f"{material_id} has {len(material.fits)} fits; fit = {fit_index} "
                         f"is out of range")
    return PartMaterial(material, fit_index)


def load_study(path, materials_dir=MATERIALS_DIR):
    path = Path(path)
    with open(path, "rb") as f:
        raw = tomllib.load(f)

    if "design" not in raw:
        raise StudyError(f"{path.name}: missing design")
    design = load_design(raw["design"])

    entries = {p.name: p.default_material for p in design.parts}
    for part, entry in raw.get("materials", {}).items():
        if part not in entries:
            raise StudyError(f"{path.name}: design {design.id} has no part {part!r}")
        entries[part] = entry
    materials = {part: _material(part, e, materials_dir) for part, e in entries.items()}
    for part, pm in materials.items():
        entry = entries[part] if isinstance(entries[part], dict) else {}
        if not pm.is_sheet:
            if "matrix" in entry or "stiffness_scale" in entry:
                raise StudyError(f"{path.name}: [materials] {part}: matrix and stiffness_scale "
                                 f"apply to sheet materials only")
            continue
        dp = design.part(part)
        if not dp.layer_thickness_mm:
            raise StudyError(f"{path.name}: sheet material {pm.material.id} on part {part!r}, "
                             f"which has no layer_thickness_mm in design {design.id}")
        if "base" not in materials and "matrix" not in entry:
            raise StudyError(f"{path.name}: [materials] {part}: give a matrix material")
        pm.matrix = (_material(f"{part}.matrix", entry["matrix"], materials_dir)
                     if "matrix" in entry else materials["base"])
        if pm.matrix.is_sheet:
            raise StudyError(f"{path.name}: [materials] {part}: the matrix must be a solid")
        pm.stiffness_scale = float(entry.get("stiffness_scale", 1.0))
        if pm.stiffness_scale <= 0:
            raise StudyError(f"{path.name}: [materials] {part}: stiffness_scale must be > 0")
        pm.layer = {"thickness_mm": float(dp.layer_thickness_mm),
                    "warp_axis": dp.warp_axis or "x", "weft_axis": dp.weft_axis or "y"}

    load = LoadCase(**raw.get("load", {}))
    if load.case not in LOAD_CASES:
        raise StudyError(f"{path.name}: load.case must be one of {LOAD_CASES}")
    if load.levels_kPa:
        lv = [float(v) for v in load.levels_kPa]
        if lv[0] <= 0 or any(b <= a for a, b in zip(lv, lv[1:])):
            raise StudyError(f"{path.name}: levels_kPa must be positive and increasing")
        load.pressure_kPa, load.steps = lv[-1], len(lv)
    elif load.pressure_kPa <= 0 or load.steps < 1:
        raise StudyError(f"{path.name}: load needs pressure_kPa > 0 and steps >= 1")

    return Study(
        name=raw.get("name", path.stem),
        path=path,
        design=design,
        materials=materials,
        mesh=MeshSpec(**raw.get("mesh", {})),
        load=load,
        backend=raw.get("backend", "febio"),
        bulk_to_shear=float(raw.get("bulk_to_shear", 100.0)),
        solver=raw.get("solver", {}),
        contact=raw.get("contact", {}),
        raw=raw,
    )
