"""FEBio 4 backend: write .feb, run the solver, read the text logs."""

from .materials import UnsupportedMaterial, febio_initial_shear_MPa, febio_material
from .results import log_stats, log_status, read_node_data, read_results
from .runner import FEBioNotFound, febio_exe, require_febio, run_febio
from .writer import LOG_FILES, write_feb

INPUT_NAME = "model.feb"


def write(run_dir, mesh, study, bend_dof, end_time=1.0):
    return write_feb(f"{run_dir}/{INPUT_NAME}", mesh, study.materials, study.load,
                     study.bulk_to_shear, study.solver, bend_dof=bend_dof,
                     contact=study.contact, end_time=end_time)


def solve(run_dir, timeout_s=None, threads=0, input_name=INPUT_NAME):
    info = run_febio(f"{run_dir}/{input_name}", timeout_s=timeout_s, threads=threads)
    info["stats"] = log_stats(info["log"])
    return info


def available():
    return febio_exe() is not None


def status(run_dir):
    return log_status(f"{run_dir}/model.log")


def read(run_dir, n_nodes):
    return read_results(run_dir, LOG_FILES, n_nodes)


__all__ = ["write", "solve", "read", "status", "available", "write_feb", "run_febio",
           "read_results", "read_node_data", "febio_material", "febio_initial_shear_MPa",
           "FEBioNotFound", "UnsupportedMaterial", "require_febio", "febio_exe"]
