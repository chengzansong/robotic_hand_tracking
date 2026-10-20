# Simulation tests

Pytest suite for the simulation module. Run it from the repository root:

```bash
.venv\Scripts\python -m pytest simulation/tests -q
```

It covers:

- loading and converting material cards
- meshing the finger CAD and tagging its surfaces
- writing FEBio input files
- parsing results and computing metrics
- run-time estimates

One test solves a single element in FEBio and compares it with the analytic solution. It is skipped when FEBio is not configured; see [../README.md](../README.md).

Cross-module tests, such as perception to control to simulation, will go in a top-level `tests/` folder once those modules exchange data.
