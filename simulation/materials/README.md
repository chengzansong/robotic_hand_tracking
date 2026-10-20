# Material cards

One TOML file per material. Studies refer to a material by its file name without `.toml`, for example `ecoflex_00-30`. Swapping a material is a one-line change in a study or a design. Both stages read these cards through [simulation/core/materials.py](../core/materials.py).

Files starting with `_` are for testing only and must never back a design decision.

## Schema

```toml
id = "ecoflex_00-30"                 # must match the file name
name = "Smooth-On Ecoflex 00-30"
family = "ecoflex"                   # ecoflex, dragon_skin, ...
density_kg_m3 = 1070

[datasheet]                          # manufacturer numbers, for reference and sanity checks
shore_hardness = "00-30"
modulus_100_kPa = 69                 # stress at 100% strain (optional)
elongation_at_break_pct = 900
tensile_strength_kPa = 1380
source = "Smooth-On Ecoflex series technical bulletin"
url = "https://..."

[purchase]                           # budget comparison; prices drift, so always date them
kit_size = "Trial unit, 2 lb (0.9 kg)"
kit_mass_kg = 0.91
price_usd = 35.00
price_per_kg_usd = 38.5
vendor = "Smooth-On store"
url = "https://..."
checked = "2026-10-19"
notes = "Also sold by Reynolds Advanced Materials and Amazon; check shipping"

[[fits]]                             # one or more constitutive fits; the first is the default
model = "ogden"                      # neo_hookean | mooney_rivlin | yeoh | ogden
mu_kPa = [ ... ]                     # ogden only
alpha = [ ... ]                      # ogden only
test = "uniaxial tension, ASTM D412"
stretch_range = [1.0, 6.0]           # range of stretch the fit was made on
source = "Marechal et al. 2021, Soft Robotics, Table III"
doi = "10.1089/soro.2019.0115"
confidence = "medium"                # high | medium | low
notes = ""
```

## Parameter conventions

Every card uses these conventions, which match Abaqus and most papers. Solver backends convert from them. FEBio's Ogden `c_i` equals `2 * mu_i`, for example. The conversion lives in the backend, never in the card.

| `model` | Keys (kPa) | Strain energy, isochoric part | Initial shear modulus μ0 |
| --- | --- | --- | --- |
| `neo_hookean` | `mu_kPa` | `μ/2 (Ī1 − 3)` | `μ` |
| `mooney_rivlin` | `c10_kPa`, `c01_kPa` | `C10 (Ī1 − 3) + C01 (Ī2 − 3)` | `2 (C10 + C01)` |
| `yeoh` | `c10_kPa`, `c20_kPa`, `c30_kPa` | `Σ Ci0 (Ī1 − 3)^i` | `2 C10` |
| `ogden` | `mu_kPa = [...]`, `alpha = [...]` | `Σ 2 μi / αi² (λ̄1^αi + λ̄2^αi + λ̄3^αi − 3)` | `Σ μi` |

Volumetric behavior is near-incompressible. A fit may give `bulk_modulus_kPa` or `poisson`. If it gives neither, the study's `bulk_to_shear` ratio applies.

In this canonical Ogden form, each term is stable when `mu_i > 0`. The often-quoted `mu_i · alpha_i > 0` rule belongs to the original Ogden form.

When a paper uses another convention, convert it into this one and say so in `notes`. Never invent coefficients. If only datasheet numbers exist, a `neo_hookean` fit estimated from `modulus_100_kPa` is allowed with `confidence = "low"`, and `notes` must show the arithmetic.

## State of the library (2026-10-19)

Eight cards: Ecoflex 00-10, 00-20, 00-30, 00-50; Dragon Skin 10, 20, 30; Elastosil M4601. All Smooth-On grades cost $33–39 per 2 lb Trial Unit. Each card has three kinds of fit, which bracket the real stiffness:

1. **Default:** the best-supported fit. For Ecoflex it is a Yeoh fit converted from Liao et al. 2020. For Dragon Skin it is our Ogden refit of the Marechal et al. 2021 raw data.
2. **Alternate:** an Ogden or Yeoh refit, labelled `low`. The Ecoflex Ogden refits are probably too soft, because the source measured strain from crosshead travel.
3. **Datasheet neo-Hookean** (`mu = P_100 / 1.75`): a common basis across grades.

Stiffness for one grade spreads 2–4× across these sources. Run Stage 1 studies across all three fits, and treat our own dogbone tests (bench plan P0.3) as the tie-breaker. Marechal et al. 2021 Table III coefficients are **not** used; each card's header explains why.

## Tools

- `tools/refit_marechal.py <Tensile-Tests-Data dir>` reproduces the refits from the Marechal raw data. The data comes from Zenodo, doi 10.5281/zenodo.3611328, and is not stored in this repo.
- `tools/validate_cards.py <repo root> <Tensile-Tests-Data dir>` checks every card: initial shear modulus, stress at λ = 2 against the datasheet, uniaxial and equibiaxial stability, and price per kg.
