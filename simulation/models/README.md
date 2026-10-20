# Physical models (finger designs)

One folder per finger design, shared by both stages. Each folder holds the source CAD and a `design.toml` that names the parts and their default materials. [simulation/core/designs.py](../core/designs.py) loads them.

```
models/
  finger_poc/         proof-of-concept PneuNet finger from Onshape
    design.toml       parts, roles, default material per part, fem_tags (surface selectors)
    cad/
      assembly.step   top + base as positioned in Onshape (the mesher uses this)
      top.step        chambered body
      base.step       base slab
  tools/
    make_pneunet.py   parametric generator for thickness variants
  pn_*/               generated variants (not in git, see below)
```

## finger_poc

The finger is 112 × 15 × 19 mm. In the assembly frame, the base occupies z = 0–4 mm and the top body z = 3–19 mm.

- **Top body:** 11 chambers, each 6 × 11 × 13 mm inside, on a 10 mm pitch. The walls between chambers are 1 mm, the open gaps between neighbouring walls are 2 mm, the side walls and roof are 2 mm, and the end walls are 3 mm. Neighbouring chambers connect through a 2 × 2 mm channel in a 3 mm bridge under each gap.
- **Base:** a flat 4 mm slab.
- **Lip and pegs:** the top has a 1 mm-wide rim and 20 small pegs on its underside that sink 1 mm into the base (338 mm³). They get pressed into the wet base silicone during casting, so the cured finger is a single piece. The FEM treats the overlap as solid top material bonded to the base.
- **Inlet:** the CAD has no inlet port. The real finger will be pierced through the side for the tube after casting.

## Parametric variants

`tools/make_pneunet.py` rebuilds this finger from its dimensions; the default parameters reproduce `finger_poc` exactly. It writes a new design folder for a variant:

```bash
# one variant
.venv\Scripts\python -m simulation.models.tools.make_pneunet pn_base3_wall1 --base_mm 3

# the standard set used in the experiment log
.venv\Scripts\python -m simulation.models.tools.make_pneunet --standard-variants
```

| Design | Change from finger_poc |
| --- | --- |
| `pn_base4_wall1` | none (exact copy, used as a reference) |
| `pn_base2_wall1`, `pn_base3_wall1`, `pn_base6_wall1` | base 2 / 3 / 6 mm |
| `pn_base4_wall0p75`, `pn_base4_wall1p5` | walls between chambers 0.75 / 1.5 mm. Pitch and gap stay fixed, and the chamber cavity absorbs the change |

Generated folders are about 1.5 MB of STEP each, so they are not stored in git. Regenerate them with the command above.

## Adding a design by hand

Copy `finger_poc/`, replace the STEP files, and change `id` in `design.toml` to match the new folder name. The expected format is STEP (AP242) from Onshape. The mesher matches the parts in the assembly to the part STEP files by volume, so the top and base must have different volumes. Meshes and solver files go to `simulation/outputs/`, never into this folder.
