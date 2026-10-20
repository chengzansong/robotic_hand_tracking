"""Regenerate renders for an existing run folder without re-solving.

    python -m simulation.stage1_finger.render simulation/outputs/<study>/<run>
    python -m simulation.stage1_finger.render <run> --center-sections

--center-sections renders only the exact centre-plane sections at every saved
pressure level, with GIFs, into <run>/renders/center/.
"""

import argparse

from .analysis.render import render_center_sections, render_run


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--center-sections", action="store_true",
                    help="only the centre-plane sections and animations (renders/center/)")
    args = ap.parse_args(argv)
    if args.center_sections:
        paths, summary = render_center_sections(args.run_dir)
        for row in summary:
            gaps = {s: g for s, g in row["slot_min_gap_on_plane_mm"].items() if g is not None}
            closest = min(gaps.items(), key=lambda kv: kv[1]) if gaps else None
            print(f"p = {row['pressure_kPa']:g} kPa: {row['contact_nodes']} on-plane wall nodes "
                  f"in contact, {row['near_nodes']} near, inverted elements {row['inverted']}"
                  + (f", closest {closest[0]} {closest[1]:.3f} mm" if closest else ""))
    else:
        paths = render_run(args.run_dir)
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
