"""Shared fixtures for the Stage 1 simulation tests. Run from the repo root:

    .venv\\Scripts\\python -m pytest simulation/tests -q
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

STUDIES = REPO_ROOT / "simulation" / "stage1_finger" / "studies"


@pytest.fixture(scope="session")
def smoke_mesh(tmp_path_factory):
    """Coarse finger_poc mesh (the smoke study's size), built once per session."""
    from simulation.core.designs import load_design
    from simulation.stage1_finger.mesh import MeshSpec, build_mesh
    out = tmp_path_factory.mktemp("mesh")
    mesh, report = build_mesh(load_design("finger_poc"), MeshSpec(size_mm=4.0, wall_size_mm=1.5),
                              out_dir=out)
    return mesh, report, out


@pytest.fixture(scope="session")
def half_mesh():
    """Coarse half model (y mid-plane symmetry)."""
    from simulation.core.designs import load_design
    from simulation.stage1_finger.mesh import MeshSpec, build_mesh
    return build_mesh(load_design("finger_poc"),
                      MeshSpec(size_mm=4.0, wall_size_mm=1.5, symmetry="half_y"))
