# Stage 1: single-finger FEM design validation

Checks a finger design before casting. The pipeline runs a quasi-static pressure sweep in FEM and reports:

- bend angle and tip path
- blocked tip force
- cavity volume versus pressure, V(p), which feeds the pump model
- peak wall stretch
- sensitivity to material and wall thickness

Playback will later turn PWM commands into a rendered finger by interpolating the sweep.

```
stage1_finger/
  study.py         # study TOML loader: design + material fit per part + mesh + load case
  run.py           # CLI: mesh -> .feb -> solve -> results, plots, renders
  render.py        # CLI: regenerate renders for an existing run folder
  summarize.py     # CLI: every run's results in one table (outputs/summary.md, summary.csv)
  mesh/            # STEP -> tagged tet10 mesh (gmsh); FEMesh container and element math
  backends/
    __init__.py    # BACKENDS: backend name -> module
    febio/         # materials.py (card -> FEBio), writer.py (.feb), runner.py, results.py (log parser)
  studies/         # design studies (eco30_sweep_10kPa, eco50_sweep_10kPa, poc_*, ...)
    compare/       # design comparisons on mesh A (relative only): baseline, base material, chamber height, blocked force
    dev/           # solver tuning, mesh convergence, and superseded studies, kept for reproducibility
  analysis/        # metrics.py (definitions), curvature.py (turn radius), sweep.py, plots.py, render.py
  playback/        # planned: PWM -> pump model -> p(t) -> interpolated, rendered finger
```

## Run

From the repository root:

```bash
.venv\Scripts\python -m simulation.stage1_finger.run simulation/stage1_finger/studies/poc_smoke.toml --no-solve
.venv\Scripts\python -m simulation.stage1_finger.run simulation/stage1_finger/studies/eco30_sweep_10kPa.toml --estimate
.venv\Scripts\python -m simulation.stage1_finger.run simulation/stage1_finger/studies/poc_pressure_sweep.toml
.venv\Scripts\python -m simulation.stage1_finger.render simulation/outputs/<study>/<run>
.venv\Scripts\python -m simulation.stage1_finger.summarize
.venv\Scripts\python -m pytest simulation/tests -q
```

**Check the time first.** Every solve starts with an estimate probe: it solves the first pressure level at full settings in `<run>/estimate/`, extrapolates the full run as a range (probe time per level x number of levels, times 1 to 2.5), and refuses to start if the upper bound exceeds `--max-minutes` (default 20). `--estimate` stops after the probe. Estimates and actual times are appended to `simulation/outputs/estimates.csv`. FEBio has no GPU path; `[solver] threads` sets OMP_NUM_THREADS.

Without FEBio configured, a run still meshes, writes `model.feb`, metadata, and `renders/mesh_tags.png`, then exits with code 2 and says how to set `febio_exe`. The FEBio uniaxial test is skipped.

Each run writes `simulation/outputs/<study>/<UTC timestamp>/` (gitignored), containing `study.toml`, `mesh/` (`mesh.msh`, `mesh.vtu`, `mesh.npz`), `model.feb`, `metadata.json` (material fits and solver parameters used, mesh stats and quality, geometry findings, git commit), and, after a solve, `model.log`, the node_data logs, `results.json`, `displacements.npz`, `plots/`, and `renders/`.

## Study format

```toml
design = "finger_poc"
backend = "febio"
bulk_to_shear = 100.0        # K / mu0 when the fit gives neither bulk_modulus_kPa nor poisson
[materials]                  # per part; omitted parts use design defaults
top = "ecoflex_00-30"                          # fit 0
base = { material = "ecoflex_00-30", fit = 2 } # fit by index
[mesh]
size_mm = 2.0
wall_size_mm = 0.8           # refined size at the cavity walls
wall_distance_mm = 2.0
symmetry = "half_y"          # optional half model on the y mid-plane (metrics are full-model)
[load]
case = "free"                # free | blocked_force
pressure_kPa = 20.0          # evenly spaced levels ...
steps = 10
# levels_kPa = [0.5, 1, 2, 3, 5]   # ... or explicit output levels
[contact]                    # slot-wall contact, on by default
enabled = true
fric_coeff = 0.0             # also: penalty, auto_penalty, two_pass, laugon, search_radius,
                             # search_tol, tolerance, penetration_tol_mm (post-solve check)
[solver]                     # optional: substeps, aggressiveness, max_retries, opt_iter,
                             # max_refs, max_ups (0 = full Newton), symmetric_stiffness,
                             # pressure_symmetric, threads, dtol, etol, rtol, lstol,
                             # domain_type, elem_type
```

## Model

- Mesh: the assembly STEP is fragmented so top and base share a conformal interface, and each part is its own volume and material. Cavity, clamp, and tip are found by the geometric selectors in `[fem_tags]` of `design.toml`. The mesh is quadratic tetrahedra (tet10), refined toward the cavity walls.
- FEBio: uncoupled hyperelastic materials. Neo-Hookean and Mooney-Rivlin cards map to `Mooney-Rivlin`, Yeoh to `Yeoh`, and Ogden to `Ogden` (c_i = 2 mu_i). The domain is the default tet10 `elastic-solid`. Cavity pressure is ramped linearly with must points at each output level. The clamp and the blocked-force tip hold are prescribed zero displacements, because FEBio logs no reactions on "zero displacement" dofs. FEBio logs a reaction as the force the body exerts on the support (checked on the single-element test), so the blocked force is `R_tip . b`, positive toward the bend side.
- Contact: the facing walls of each of the 10 slots between chambers form one `sliding-elastic` pair each (two-pass, auto penalty, frictionless by default). Without it the walls of neighbouring chambers pass through each other from about 1.5 kPa with the placeholder material. After the solve, `analysis/contact.py` measures the minimum signed gap per slot and per step (results.json `slot_min_gap_mm`, warning when below `-penetration_tol_mm`) and renders `contact_gaps_max.png`.
- Not modelled: the tip curling back onto the base or clamp. `results.json` flags bends past 300 degrees or a tip within 1 mm of the proximal half (`curl`).
- Partial runs: if FEBio stops converging, the converged levels are still read, plotted, and rendered; `results.json` has `reached_full_load = false` and `last_converged_pressure_kPa`.
- Metrics: see the docstring of `analysis/metrics.py`. The bend angle is the in-plane rotation of the tip-face normal, positive toward the base.
- Turn radius (free bending only): see `analysis/curvature.py`. A circle fit to the deformed base underside on the centre plane, over the chambered span, plus circularity, mean radius (arc length / bend), the curvature profile along the finger, and the underside length change. In `results.json` as `turn_radius`, plotted to `plots/turn_radius_sides.png` (dimensioned side views with reference objects) and `plots/turn_radius.png`.

## Solver tuning (2026-10-20, smoke mesh, placeholder material)

Full model, 42k nodes / 127k equations, 2 kPa in 4 levels, no contact, one change at a time from the first solve. Bend and dV agree to 4 significant figures in every row except T6.

| Run | Change | Wall (s) | Steps | Iterations | Reforms | Bend (deg) | dV (mm^3) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| baseline | BFGS, symmetric, 4 substeps | 507 | 42 | 1249 | 169 | 51.114 | 2305.9 |
| T1 | non-symmetric solver and pressure load | 522 | 41 | 1062 | 144 | 51.114 | 2305.9 |
| T2 | full Newton (max_ups = 0) | 393* | 16 | 243 | 243 | 51.114 | 2305.9 |
| T1+T2 | both | 440 | 16 | 243 | 243 | 51.114 | 2305.9 |
| T3a | Newton, 1 substep | 303 | 10 | 171 | 171 | 51.114 | 2305.9 |
| T3b | Newton, 2 substeps, aggressiveness 1 | 258 | 11 | 176 | 176 | 51.114 | 2305.9 |
| T4 | Newton, OMP_NUM_THREADS = 12 | 338 | 16 | 243 | 243 | 51.114 | 2305.9 |
| T4 | Newton, OMP_NUM_THREADS = 24 | 325 | 16 | 243 | 243 | 51.114 | 2305.9 |
| T6 | Newton, TET10G4 integration | 327 | 16 | 243 | 243 | 51.143 | 2312.8 |
| T5 | Newton, half model (75.8k equations) | 210 | 15 | 232 | 232 | 51.006 | 2317.1 |

\* T2 shared the CPU with a pytest run; T4 at 24 threads is the clean full-Newton reference.

With contact (half model, 0.5 to 3 kPa, Newton, 2 substeps, aggressiveness 1): the symmetric tangent stalled at 2.35 kPa on "max nr of reformations" (760 s); the non-symmetric tangent finished in 341 s with no retries. Contact off on the same model took 321 s.

Chosen defaults (writer `SOLVER_DEFAULTS`): full Newton, non-symmetric, 2 substeps, aggressiveness 1, 24 threads, default TET10G8. On the smoke mesh that is about 2x faster than the first solve without contact, and the half model takes another 35 % off. Bend changes by 0.2 % and dV by 0.5 % with the half model. TET10G4 saves nothing because the linear solver is 75 % of the time.

Code here may import `simulation.core`. Nothing in `core/` or `stage2_hand/` imports from here.
