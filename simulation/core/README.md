# Core (shared by both stages)

Solver-independent code that Stage 1 and Stage 2 both import. Nothing here may import from `stage1_finger/` or `stage2_hand/`.

| File | Purpose |
| --- | --- |
| `materials.py` | Load material cards from [../materials/](../materials/) |
| `designs.py` | Load finger designs from [../models/](../models/) |
| `local.py` | Machine-specific paths (FEBio executable) from `simulation/local.toml` |
| `pump.py` | Planned: pump model, PWM duty to chamber pressure |
| `state.py` | Planned: finger state passed out of the module |

## Units

Configuration files carry explicit unit suffixes (`_kPa`, `_mm`, `_kg_m3`). Inside solvers, use **mm, N, MPa, s, tonne/mm³**, which is a consistent set, so 1 kPa = 0.001 MPa. Convert at the boundary of each backend.

CAD from Onshape arrives in metres (`cad_units = "m"` in `design.toml`). The mesher scales it to mm.
