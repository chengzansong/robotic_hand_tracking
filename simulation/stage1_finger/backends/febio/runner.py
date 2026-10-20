"""Run the FEBio solver as a subprocess and capture its output.

The executable path comes from simulation.core.local.local_setting("febio_exe"),
i.e. simulation/local.toml or the FEBIO_EXE environment variable.
Command line (FEBio user manual 4.7, section 2.3; FEBioLib/cmdoptions.cpp):
    febio4 -i model.feb -silent
-silent suppresses screen output but FEBio still writes model.log and
model.xplt next to the input file.
"""

import os
import subprocess
import time
from pathlib import Path

from simulation.core.local import local_setting

INSTALL_HINT = (
    "FEBio is not configured. Install FEBio Studio from https://febio.org/downloads/ "
    "(free registration), then set febio_exe in simulation/local.toml "
    "(copy simulation/local.example.toml) or set the FEBIO_EXE environment variable "
    "to the bundled solver, for example C:\\Program Files\\FEBioStudio\\bin\\febio4.exe."
)


class FEBioNotFound(RuntimeError):
    pass


def febio_exe():
    """Path to the FEBio executable, or None when unset or missing."""
    exe = local_setting("febio_exe")
    if not exe:
        return None
    exe = Path(exe)
    return exe if exe.exists() else None


def require_febio():
    exe = febio_exe()
    if exe is None:
        configured = local_setting("febio_exe")
        detail = f" (configured path {configured!r} does not exist)" if configured else ""
        raise FEBioNotFound(INSTALL_HINT + detail)
    return exe


def run_febio(feb_path, timeout_s=None, threads=0):
    """Solve feb_path in its own folder. Returns a dict with paths and return code.

    threads > 0 sets OMP_NUM_THREADS for the solver (FEBio and its Pardiso
    solver use OpenMP; FEBio does not print the thread count it uses).
    A timeout kills the solver and returns returncode "timeout".
    """
    exe = require_febio()
    feb_path = Path(feb_path).resolve()
    stdout_path = feb_path.with_name("febio_stdout.txt")
    env = dict(os.environ)
    if threads:
        env["OMP_NUM_THREADS"] = str(threads)
    t0 = time.time()
    with open(stdout_path, "w", encoding="utf-8", errors="replace") as out:
        try:
            rc = subprocess.run([str(exe), "-i", feb_path.name, "-silent"],
                                cwd=feb_path.parent, stdout=out, stderr=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL, timeout=timeout_s, env=env).returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
    return {
        "exe": str(exe),
        "threads": threads or "unset",
        "returncode": rc,
        "wall_time_s": round(time.time() - t0, 2),
        "log": str(feb_path.with_suffix(".log")),
        "plot": str(feb_path.with_suffix(".xplt")),
        "stdout": str(stdout_path),
    }
