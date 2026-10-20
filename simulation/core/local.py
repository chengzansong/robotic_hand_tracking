"""Machine-specific settings, such as solver executable paths.

Read from simulation/local.toml, which is gitignored. Copy
simulation/local.example.toml to start. An environment variable of the
same name in upper case (for example FEBIO_EXE) takes precedence.
"""

import os
import tomllib
from pathlib import Path

LOCAL_FILE = Path(__file__).resolve().parent.parent / "local.toml"


def local_setting(key, default=None):
    env = os.environ.get(key.upper())
    if env:
        return env
    if LOCAL_FILE.exists():
        with open(LOCAL_FILE, "rb") as f:
            return tomllib.load(f).get(key, default)
    return default
