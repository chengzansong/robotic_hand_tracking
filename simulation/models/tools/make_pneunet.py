"""Generate a parametric PneuNet finger design folder (STEP parts + design.toml).

The defaults reproduce finger_poc exactly (top 12340 mm^3, base 6720 mm^3,
lip and pegs 338 mm^3, cavity 9598 mm^3), measured from the Onshape STEP by
cross-section on 2026-10-19. Variants change one dimension at a time:

- base_mm: base slab thickness
- wall_mm: wall between neighbouring chambers. Pitch and slot stay fixed, so
  finger length and the contact gap are unchanged; the chamber cavity length
  absorbs the change (cavity_x = pitch - slot - 2 * wall).
- side_mm, roof_mm: side walls and roof over the chambers

Assembly frame (mm): finger axis +x from 0, width along y from 0, base on
z = 0..base_mm, chambered top above it with a 1 mm lip and pegs sunk into the
base, the same convention as finger_poc.

Usage, from the repository root:
    python -m simulation.models.tools.make_pneunet finger_b3 --base_mm 3
    python -m simulation.models.tools.make_pneunet --standard-variants
"""

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

import gmsh

MODELS_DIR = Path(__file__).resolve().parent.parent


@dataclass
class PneuNet:
    n_chambers: int = 11
    pitch_mm: float = 10.0
    slot_mm: float = 2.0          # open gap between neighbouring chamber walls
    wall_mm: float = 1.0          # wall between neighbouring chambers
    end_mm: float = 3.0           # proximal and distal end walls
    width_mm: float = 15.0
    side_mm: float = 2.0
    cavity_h_mm: float = 13.0     # chamber height above the base top
    roof_mm: float = 2.0
    bridge_h_mm: float = 3.0      # solid under each slot, holding the feed channel
    channel_w_mm: float = 2.0
    channel_h_mm: float = 2.0
    base_mm: float = 4.0
    lip_mm: float = 1.0           # lip and pegs sunk into the base
    lip_w_mm: float = 1.0         # rim width, under the inner edge of the end and side walls
    peg_y_mm: tuple = (2.0, 4.5)  # peg extent from the centreline, both sides

    @property
    def cavity_x_mm(self):
        return self.pitch_mm - self.slot_mm - 2 * self.wall_mm

    @property
    def length_mm(self):
        # n cavities, n-1 (wall, slot, wall) gaps, two end walls
        return 2 * self.end_mm + self.n_chambers * self.pitch_mm - self.slot_mm - 2 * self.wall_mm

    def check(self):
        if self.cavity_x_mm <= 0:
            raise ValueError("walls too thick for the pitch and slot")
        if self.lip_mm >= self.base_mm:
            raise ValueError("lip must be thinner than the base")


def _box(x0, y0, z0, x1, y1, z1):
    return gmsh.model.occ.addBox(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0)


def _chamber_x(p):
    """(x0, x1) of each chamber cavity."""
    out, x = [], p.end_mm
    for _ in range(p.n_chambers):
        out.append((x, x + p.cavity_x_mm))
        x += p.cavity_x_mm + 2 * p.wall_mm + p.slot_mm
    return out


def _top(p):
    occ = gmsh.model.occ
    z0, w, L = p.base_mm, p.width_mm, p.length_mm
    zt = z0 + p.cavity_h_mm + p.roof_mm
    body = _box(0, 0, z0, L, w, zt)
    cuts = []
    chambers = _chamber_x(p)
    for x0, x1 in chambers:
        cuts.append(_box(x0, p.side_mm, z0, x1, w - p.side_mm, z0 + p.cavity_h_mm))
    yc = w / 2
    for (a0, a1), (b0, b1) in zip(chambers, chambers[1:]):
        s0, s1 = a1 + p.wall_mm, b0 - p.wall_mm
        cuts.append(_box(s0, 0, z0 + p.bridge_h_mm, s1, w, zt))
        cuts.append(_box(a1, yc - p.channel_w_mm / 2, z0, b0, yc + p.channel_w_mm / 2,
                         z0 + p.channel_h_mm))
    body = occ.cut([(3, body)], [(3, c) for c in cuts])[0]
    # lip: a rim under the side and end walls, plus two pegs under each bridge
    zl = z0 - p.lip_mm
    ex, sy, wl = p.end_mm, p.side_mm, p.lip_w_mm
    rim_out = _box(ex - wl, sy - wl, zl, L - ex + wl, w - sy + wl, z0)
    rim_in = _box(ex, sy, zl, L - ex, w - sy, z0)
    rim = occ.cut([(3, rim_out)], [(3, rim_in)])[0]
    pegs = []
    for (a0, a1), (b0, b1) in zip(chambers, chambers[1:]):
        s0, s1 = a1 + p.wall_mm, b0 - p.wall_mm
        for lo, hi in ((yc + p.peg_y_mm[0], yc + p.peg_y_mm[1]),
                       (yc - p.peg_y_mm[1], yc - p.peg_y_mm[0])):
            pegs.append((3, _box(s0, lo, zl, s1, hi, z0)))
    top = occ.fuse(body, rim + pegs)[0]
    occ.synchronize()
    return top


def _base(p):
    t = _box(0, 0, 0, p.length_mm, p.width_mm, p.base_mm)
    gmsh.model.occ.synchronize()
    return [(3, t)]


def _write(builder, p, path):
    gmsh.clear()
    ents = builder(p)
    gmsh.write(str(path))
    return sum(gmsh.model.occ.getMass(d, t) for d, t in ents)


def _write_assembly(p, path):
    gmsh.clear()
    _top(p)
    _base(p)
    gmsh.write(str(path))


def make(design_id, p, material="ecoflex_00-30", note=""):
    p.check()
    out = MODELS_DIR / design_id
    (out / "cad").mkdir(parents=True, exist_ok=True)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        v_top = _write(_top, p, out / "cad" / "top.step")
        v_base = _write(_base, p, out / "cad" / "base.step")
        _write_assembly(p, out / "cad" / "assembly.step")
    finally:
        gmsh.finalize()
    if abs(v_top - v_base) < 1e-3 * v_base:
        raise ValueError("top and base volumes are equal; the mesher matches parts by volume")

    cavity = p.n_chambers * p.cavity_x_mm * (p.width_mm - 2 * p.side_mm) * p.cavity_h_mm \
        + (p.n_chambers - 1) * (p.slot_mm + 2 * p.wall_mm) * p.channel_w_mm * p.channel_h_mm
    params = "\n".join(f"{k} = {list(v) if isinstance(v, tuple) else v}"
                       for k, v in asdict(p).items())
    (out / "design.toml").write_text(f'''# Generated by simulation/models/tools/make_pneunet.py. Edit the generator, not this file.
id = "{design_id}"
name = "Parametric PneuNet {design_id}"
cad_units = "mm"
assembly = "cad/assembly.step"
notes = """{note}
Generated volumes (mm^3): top {v_top:.1f}, base {v_base:.1f}, cavity {cavity:.1f}.
Length {p.length_mm:.1f} mm, chamber cavity length {p.cavity_x_mm:.2f} mm."""

[[parts]]
name = "top"
step = "cad/top.step"
role = "chambered body"
default_material = "{material}"

[[parts]]
name = "base"
step = "cad/base.step"
role = "base slab, strain limiting by thickness"
default_material = "{material}"

[fem_tags]
axis = "x"
clamp_end = "min"
bend_toward = "-z"
plane_tol_mm = 0.01
cavity = "inner_shells"
interior_margin_mm = 0.5
overlap_owner = "top"

[parametric]
{params}
''', encoding="utf-8")
    return {"top_mm3": v_top, "base_mm3": v_base, "cavity_mm3": cavity, "length_mm": p.length_mm}


# The thickness variants used in the experiment log. Generated folders are
# gitignored; rebuild them with --standard-variants.
STANDARD_VARIANTS = {
    "pn_base4_wall1": {},
    "pn_base2_wall1": {"base_mm": 2.0},
    "pn_base3_wall1": {"base_mm": 3.0},
    "pn_base6_wall1": {"base_mm": 6.0},
    "pn_base4_wall0p75": {"wall_mm": 0.75},
    "pn_base4_wall1p5": {"wall_mm": 1.5},
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("design_id", nargs="?")
    for f in ("base_mm", "wall_mm", "side_mm", "roof_mm"):
        ap.add_argument(f"--{f}", type=float)
    ap.add_argument("--material", default="ecoflex_00-30")
    ap.add_argument("--standard-variants", action="store_true",
                    help="generate every design in STANDARD_VARIANTS")
    a = ap.parse_args()
    if a.standard_variants:
        for design_id, dims in STANDARD_VARIANTS.items():
            print(design_id, make(design_id, PneuNet(**dims), a.material))
        return
    if not a.design_id:
        ap.error("give a design_id or --standard-variants")
    p = PneuNet(**{f: getattr(a, f) for f in ("base_mm", "wall_mm", "side_mm", "roof_mm")
                   if getattr(a, f) is not None})
    print(make(a.design_id, p, a.material))


if __name__ == "__main__":
    main()
