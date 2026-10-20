# Simulation

Takes pump commands from [control/](../control/) and a physical model of the finger, and outputs a rendered simulation plus the predicted finger or hand state. See [CONTRACTS.md](../CONTRACTS.md).

The work is split into two stages. They share data and core code, but neither imports the other.

| Stage | Folder | Job | Tools | State |
| --- | --- | --- | --- | --- |
| 1 | [stage1_finger/](stage1_finger/) | One finger, before casting: FEM of bend angle, tip force, cavity volume, and wall stretch against pressure | gmsh, FEBio, PyVista | Working. Results in [docs/experiments.md](docs/experiments.md) |
| 2 | [stage2_hand/](stage2_hand/) | Whole hand, control, and contact with objects | MuJoCo, calibrated from Stage 1 | Not started |

## Layout

```
simulation/
  core/            shared, solver independent: material and design loaders, local paths
  materials/       one TOML card per material, with fits, sources, and prices (see materials/README.md)
  models/          one folder per finger design: design.toml + cad/*.step
    finger_poc/    the proof-of-concept finger from Onshape
    tools/         make_pneunet.py: parametric finger generator (thickness variants)
  stage1_finger/   Stage 1: meshing, FEBio backend, studies, analysis, renders
  stage2_hand/     Stage 2 (not started)
  docs/            experiment log and key images
  tests/           pytest suite (run: python -m pytest simulation/tests -q)
  outputs/         generated runs (not in git)
  local.toml       machine paths such as the FEBio executable (not in git)
```

Swapping parts happens in configuration files, not code. A **study** (a TOML file in `stage1_finger/studies/`) names a design, a material and fit for each part, the mesh size, the pressure levels, contact, and solver settings.

| To add | Do this |
| --- | --- |
| A material | Add a card to `materials/` |
| A finger design | Add a folder to `models/`, or generate one with `models/tools/make_pneunet.py` |
| A different FEM solver | Add it under `stage1_finger/backends/` |

## Setup

From the repository root, with Python 3.11:

```bash
py -3.11 -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
copy simulation\local.example.toml simulation\local.toml
```

FEBio is not a pip package:

1. Register for a free account at [febio.org](https://febio.org/downloads/) and install **FEBio Studio**. The SDK is not needed.
2. Find the bundled solver `febio4.exe`; this is not the `FEBioStudio.exe` GUI. The default install puts it at `C:\Program Files\FEBioStudio\bin\febio4.exe`.
3. Put that path in `simulation/local.toml` as `febio_exe = '...'` (single quotes), or set the `FEBIO_EXE` environment variable.

Stage 1 needs no GPU. FEBio runs on the CPU only.

To check the install:

```bash
.venv\Scripts\python -m pytest simulation/tests -q
```

This includes a one-element FEBio check against the analytic solution.

## Running a study

```bash
# 1. Estimate run time first: meshes, then solves the first pressure level and extrapolates
.venv\Scripts\python -m simulation.stage1_finger.run simulation/stage1_finger/studies/poc_smoke.toml --estimate

# 2. Run it, with a hard timeout in seconds
.venv\Scripts\python -m simulation.stage1_finger.run simulation/stage1_finger/studies/poc_smoke.toml --timeout-s 1800

# 3. Look at the renders, and redraw them later without re-solving if needed
.venv\Scripts\python -m simulation.stage1_finger.render simulation/outputs/poc_smoke/<run folder>
```

Each run writes to `simulation/outputs/<study>/<UTC timestamp>/`:

- `results.json`: bend angle, tip position, cavity volume, stretch, wall gaps, and reactions for each pressure level
- `metadata.json`: the exact material fits, mesh, and solver settings used, plus the git commit
- `plots/`
- `renders/`: the deformed shape at several pressures, a cross-section, the wall-gap map, and the tagged mesh

**Always look at the renders.** Metrics alone don't show problems such as walls passing through each other.

Parametric thickness variants (`pn_*`) are generated, not stored in git. Create them before running studies that use them:

```bash
.venv\Scripts\python -m simulation.models.tools.make_pneunet --standard-variants
```

## Units

Configuration files use explicit unit suffixes (`_kPa`, `_mm`). Solvers work in mm, N, MPa, and s. CAD from Onshape is exported in metres and imported in mm.

## Status and results

See [docs/experiments.md](docs/experiments.md). It covers the install check, contact, solver tuning, mesh convergence, the Ecoflex 00-30 pressure sweep, a first Ecoflex 00-50 comparison, and the planned next experiments.
