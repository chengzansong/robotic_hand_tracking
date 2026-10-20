"""Regenerate renders for an existing run folder without re-solving.

    python -m simulation.stage1_finger.render simulation/outputs/<study>/<run>
"""

import argparse

from .analysis.render import render_run


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir")
    args = ap.parse_args(argv)
    for path in render_run(args.run_dir):
        print(path)


if __name__ == "__main__":
    main()
