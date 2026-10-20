"""Load finger designs from simulation/models/<design>/design.toml.

A design names its CAD parts and the default material for each part.
Studies may override any part's material without touching the design.
"""

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

MODELS_DIR = Path(__file__).resolve().parent.parent / "models"


@dataclass
class Part:
    name: str
    step: Path
    role: str
    default_material: str


@dataclass
class Design:
    id: str
    name: str
    directory: Path
    cad_units: str
    parts: list[Part]
    assembly: Path | None = None
    fem_tags: dict = field(default_factory=dict)
    notes: str = ""

    def part(self, name):
        for p in self.parts:
            if p.name == name:
                return p
        raise KeyError(f"design {self.id} has no part {name!r}")

    def materials(self, overrides=None):
        """Part name -> material id, with study overrides applied."""
        chosen = {p.name: p.default_material for p in self.parts}
        for name, material_id in (overrides or {}).items():
            if name not in chosen:
                raise KeyError(f"design {self.id} has no part {name!r}")
            chosen[name] = material_id
        return chosen


def load_design(design_id, models_dir=MODELS_DIR):
    directory = Path(models_dir) / design_id
    with open(directory / "design.toml", "rb") as f:
        raw = tomllib.load(f)

    if raw.get("id") != design_id:
        raise ValueError(f"{directory.name}/design.toml: id does not match the folder name")

    parts = [
        Part(
            name=p["name"],
            step=directory / p["step"],
            role=p.get("role", ""),
            default_material=p["default_material"],
        )
        for p in raw["parts"]
    ]
    for p in parts:
        if not p.step.exists():
            raise FileNotFoundError(f"design {design_id}: missing {p.step}")

    assembly = raw.get("assembly")
    return Design(
        id=raw["id"],
        name=raw.get("name", raw["id"]),
        directory=directory,
        cad_units=raw.get("cad_units", "mm"),
        parts=parts,
        assembly=directory / assembly if assembly else None,
        fem_tags=raw.get("fem_tags", {}),
        notes=raw.get("notes", ""),
    )


def list_designs(models_dir=MODELS_DIR):
    return sorted(p.parent.name for p in Path(models_dir).glob("*/design.toml"))
