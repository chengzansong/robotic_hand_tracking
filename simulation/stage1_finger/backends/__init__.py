"""Solver backends. A study's `backend` key picks one by name.

Each backend module provides:
    write(run_dir, mesh, study, bend_dof) -> dict    solver input, written to run_dir
    available() -> bool                              solver installed and configured
    solve(run_dir, timeout_s=None) -> dict           run the solver in run_dir
    status(run_dir) -> str                           'normal', 'error', ...
    read(run_dir, n_nodes) -> (times, u, steps)      t in [0, 1], u (steps, N, 3) mm,
                                                     per-step reaction sums in N
"""

import importlib

BACKENDS = {
    "febio": "simulation.stage1_finger.backends.febio",
    # "sofa": planned alternate backend
}


def get_backend(name):
    if name not in BACKENDS:
        raise KeyError(f"unknown backend {name!r}; available: {sorted(BACKENDS)}")
    return importlib.import_module(BACKENDS[name])
