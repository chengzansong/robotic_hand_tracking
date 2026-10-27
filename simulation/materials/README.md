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

Eight cards: Ecoflex 00-10, 00-20, 00-30, 00-50; Dragon Skin 10, 20, 30; Elastosil M4601. All Smooth-On grades cost $33–39 per 2 lb Trial Unit. Most cards have three kinds of fit, which bracket the real stiffness:

1. **Default (fit 0):** the best-supported fit. For Ecoflex 00-10, 00-30 and 00-50 it is a Yeoh fit converted from Liao et al. 2020. For Ecoflex 00-20 it is the same Shore-law model, predicted rather than measured by Liao et al., so it is `low`. For Dragon Skin it is our Ogden refit of the Marechal et al. 2021 raw data. For Elastosil M4601 it is the Yeoh fit of Wang et al. 2018, the only fit on that card.
2. **Alternate:** an Ogden or Yeoh refit, labelled `low` (Ecoflex 00-10, 00-30 and 00-50 have an Ogden refit; Dragon Skin has a Yeoh refit). The Ecoflex Ogden refits are probably too soft, because the source measured strain from crosshead travel. Ecoflex 00-20 and Elastosil M4601 have no alternate.
3. **Datasheet neo-Hookean** (`mu = P_100 / 1.75`): a common basis across grades. Every card has one except Elastosil M4601.

Stiffness for one grade spreads 2–4× across these sources. Run Stage 1 studies across the fits of a card, and treat our own dogbone tests (bench plan P0.3) as the tie-breaker. Marechal et al. 2021 Table III coefficients are **not** used; each card's header explains why.

### Fits at a glance

mu0 is the initial shear modulus the fit implies (Yeoh 2 C10, Ogden sum of mu_i, neo-Hookean mu), from `core/materials.py`. Fit 0 is the default. Results reported so far used fit 0 of `ecoflex_00-30`, `ecoflex_00-50` and `dragon_skin_20` (see [Material profiles used](../docs/experiments.md#material-profiles-used)).

| Card | Fit | Model | mu0 (kPa) | Confidence | Source |
| --- | --- | --- | --- | --- | --- |
| [dragon_skin_10](dragon_skin_10.toml) | 0 | Ogden | 44.4 | medium | Library refit of Marechal et al. 2021 raw data (doi 10.1089/soro.2019.0115) |
| [dragon_skin_10](dragon_skin_10.toml) | 1 | Yeoh | 52.7 | low | Library refit of Marechal et al. 2021 raw data (doi 10.1089/soro.2019.0115) |
| [dragon_skin_10](dragon_skin_10.toml) | 2 | neo-Hookean | 86.7 | low | Smooth-On datasheet (100% modulus) |
| [dragon_skin_20](dragon_skin_20.toml) | 0 | Ogden | 125 | medium | Library refit of Marechal et al. 2021 raw data (doi 10.1089/soro.2019.0115) |
| [dragon_skin_20](dragon_skin_20.toml) | 1 | Yeoh | 148 | low | Library refit of Marechal et al. 2021 raw data (doi 10.1089/soro.2019.0115) |
| [dragon_skin_20](dragon_skin_20.toml) | 2 | neo-Hookean | 193 | low | Smooth-On datasheet (100% modulus) |
| [dragon_skin_30](dragon_skin_30.toml) | 0 | Ogden | 184 | medium | Library refit of Marechal et al. 2021 raw data (doi 10.1089/soro.2019.0115) |
| [dragon_skin_30](dragon_skin_30.toml) | 1 | Yeoh | 212 | low | Library refit of Marechal et al. 2021 raw data (doi 10.1089/soro.2019.0115) |
| [dragon_skin_30](dragon_skin_30.toml) | 2 | neo-Hookean | 339 | low | Smooth-On datasheet (100% modulus) |
| [ecoflex_00-10](ecoflex_00-10.toml) | 0 | Yeoh | 11.2 | low | Liao, Hossain, Yao 2020, doi 10.1016/j.mechmat.2020.103366 |
| [ecoflex_00-10](ecoflex_00-10.toml) | 1 | Ogden | 6.17 | low | Library refit of Marechal et al. 2021 raw data (doi 10.1089/soro.2019.0115) |
| [ecoflex_00-10](ecoflex_00-10.toml) | 2 | neo-Hookean | 31.5 | low | Smooth-On datasheet (100% modulus) |
| [ecoflex_00-20](ecoflex_00-20.toml) | 0 | Yeoh | 18.2 | low | Liao, Hossain, Yao 2020, doi 10.1016/j.mechmat.2020.103366 (Shore-law prediction) |
| [ecoflex_00-20](ecoflex_00-20.toml) | 1 | neo-Hookean | 31.5 | low | Smooth-On datasheet (100% modulus) |
| [ecoflex_00-30](ecoflex_00-30.toml) | 0 | Yeoh | 25.2 | medium | Liao, Hossain, Yao 2020, doi 10.1016/j.mechmat.2020.103366 |
| [ecoflex_00-30](ecoflex_00-30.toml) | 1 | Ogden | 14.6 | low | Library refit of Marechal et al. 2021 raw data (doi 10.1089/soro.2019.0115) |
| [ecoflex_00-30](ecoflex_00-30.toml) | 2 | neo-Hookean | 39.4 | low | Smooth-On datasheet (100% modulus) |
| [ecoflex_00-50](ecoflex_00-50.toml) | 0 | Yeoh | 39.1 | medium | Liao, Hossain, Yao 2020, doi 10.1016/j.mechmat.2020.103366 |
| [ecoflex_00-50](ecoflex_00-50.toml) | 1 | Ogden | 19.5 | low | Library refit of Marechal et al. 2021 raw data (doi 10.1089/soro.2019.0115) |
| [ecoflex_00-50](ecoflex_00-50.toml) | 2 | neo-Hookean | 47.3 | low | Smooth-On datasheet (100% modulus) |
| [elastosil_m4601](elastosil_m4601.toml) | 0 | Yeoh | 250 | low | Wang, Ge, Gu 2018, doi 10.1016/j.sna.2018.01.018 |

## Tools

- `tools/refit_marechal.py <Tensile-Tests-Data dir>` reproduces the refits from the Marechal raw data. The data comes from Zenodo, doi 10.5281/zenodo.3611328, and is not stored in this repo.
- `tools/validate_cards.py <repo root> <Tensile-Tests-Data dir>` checks every card: initial shear modulus, stress at λ = 2 against the datasheet, uniaxial and equibiaxial stability, and price per kg.
