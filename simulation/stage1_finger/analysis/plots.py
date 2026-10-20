"""Pressure-sweep plots from results.json rows (matplotlib, no display needed)."""

from pathlib import Path


def plot_sweep(rows, out_dir, title=""):
    """Write p_vs_angle.png and sweep.png (angle, V(p), stretch, force). Returns paths."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    p = [r["pressure_kPa"] for r in rows]
    angle = [r["bend_angle_deg"] for r in rows]

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(p, angle, "o-")
    ax.set(xlabel="chamber pressure (kPa)", ylabel="bend angle (deg)", title=title)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    single = out_dir / "p_vs_angle.png"
    fig.savefig(single, dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    axes[0, 0].plot(p, angle, "o-")
    axes[0, 0].set(ylabel="bend angle (deg)")
    axes[0, 1].plot(p, [r["cavity_volume_mm3"] for r in rows], "o-")
    axes[0, 1].set(ylabel="cavity volume (mm$^3$)")
    for part in rows[0]["principal_stretch"]:
        axes[1, 0].plot(p, [r["principal_stretch"][part]["max"] for r in rows], "o-", label=part)
    axes[1, 0].set(ylabel="peak principal stretch")
    axes[1, 0].legend()
    if "blocked_force_N" in rows[0]:
        axes[1, 1].plot(p, [r["blocked_force_N"] for r in rows], "o-")
        axes[1, 1].set(ylabel="blocked tip force (N)")
    else:
        tip = [r["tip_displacement_mm"] for r in rows]
        axes[1, 1].plot([t[0] for t in tip], [t[2] for t in tip], "o-")
        axes[1, 1].set(xlabel="tip dx (mm)", ylabel="tip dz (mm)", title="tip path")
        axes[1, 1].set_aspect("equal", adjustable="datalim")
    for a in axes.flat:
        a.grid(alpha=0.3)
        if not a.get_xlabel():
            a.set_xlabel("chamber pressure (kPa)")
    fig.suptitle(title)
    fig.tight_layout()
    combined = out_dir / "sweep.png"
    fig.savefig(combined, dpi=150)
    plt.close(fig)
    return [single, combined]
