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


# Target objects drawn for scale on the turn-radius side views (see
# plans/2026-10-19_finger-targets.md): (label, kind, size_mm)
REFERENCE_OBJECTS = (("grape 20", "circle", 20.0), ("cube 26", "square", 26.0),
                     ("tennis ball 67", "circle", 67.0))


def plot_turn_radius(mesh, u_all, pressures_kPa, rows, out_dir, title="",
                     axis="x", bend_toward="-z"):
    """Dimensioned side views with the fitted turn circle, and R vs pressure.

    Writes turn_radius_sides.png (one mm-scaled panel per level: the deformed
    centre-plane section, the base underside, the fitted circle with its
    radius, and reference objects at the circle centre) and turn_radius.png
    (R_fit, R_avg, circularity, curvature profile). Returns paths.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Circle, Rectangle

    from ..mesh.fe_mesh import tri6_to_linear
    from .curvature import underside_centerline
    from .metrics import axis_vector

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    a, b = axis_vector(axis), axis_vector(bend_toward)
    ids = underside_centerline(mesh, axis, bend_toward)
    by_p = {round(r["pressure_kPa"], 6): r for r in rows}
    levels = [(p, u) for p, u in zip(pressures_kPa, u_all) if round(float(p), 6) in by_p]
    sym = mesh.surfaces.get("symmetry")

    n = len(levels)
    cols = min(3, n)
    nrows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(nrows, cols, figsize=(5.2 * cols, 4.6 * nrows), squeeze=False)
    for k, (p, u) in enumerate(levels):
        ax = axes.flat[k]
        r = by_p[round(float(p), 6)]
        x = mesh.nodes + u[: len(mesh.nodes)]
        plane = lambda pts: np.column_stack([pts @ a, -(pts @ b)])   # height up, palm down
        if sym is not None:
            tris = tri6_to_linear(sym)
            ax.add_collection(PolyCollection(plane(x[tris].reshape(-1, 3)).reshape(-1, 3, 2),
                                             facecolor="0.82", edgecolor="0.7", lw=0.1))
        line = plane(x[ids])
        ax.plot(line[:, 0], line[:, 1], "-", color="tab:blue", lw=2, label="base underside")
        c = np.array([r["center_mm"][0], -r["center_mm"][1]])
        R = r["R_fit_mm"]
        style = "-" if r["fit_reliable"] else ":"
        ax.add_patch(Circle(c, R, fill=False, ls=style, color="tab:red", lw=1.2))
        ax.plot(*c, "+", color="tab:red")
        ax.plot([c[0], c[0] + R], [c[1], c[1]], color="tab:red", lw=0.8)
        ax.annotate(f"R = {R:.1f} mm", (c[0] + R / 2, c[1]), textcoords="offset points",
                    xytext=(0, 4), ha="center", color="tab:red", fontsize=9)
        if r["fit_reliable"]:
            for (label, kind, size), col in zip(REFERENCE_OBJECTS, ("tab:purple", "tab:green", "tab:olive")):
                if kind == "circle":
                    ax.add_patch(Circle(c, size / 2, fill=False, ls="--", color=col, lw=0.9, label=label))
                else:
                    ax.add_patch(Rectangle(c - size / 2, size, size, fill=False, ls="--",
                                           color=col, lw=0.9, label=label))
        pts = np.vstack([plane(x), c + [[-R, -R], [R, R]] if r["fit_reliable"] else plane(x)])
        lo, hi = pts.min(0) - 5, pts.max(0) + 5
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_aspect("equal")
        ax.grid(alpha=0.3)
        ax.set(xlabel="x (mm)", ylabel="z (mm)",
               title=f"{p:g} kPa: bend {r['span_bend_deg']:.0f}° over span, "
                     f"rms {r['circularity_rms_pct']:.2f}% of R"
                     + ("" if r["fit_reliable"] else " (fit unreliable)"))
        ax.title.set_fontsize(9)
        if k == n - 1:
            ax.legend(fontsize=7, loc="best")
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.suptitle(f"{title}: turn radius of the base underside (centre plane)", fontsize=11)
    fig.tight_layout()
    sides = out_dir / "turn_radius_sides.png"
    fig.savefig(sides, dpi=140)
    plt.close(fig)

    p = [r["pressure_kPa"] for r in rows]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    ok = [r["fit_reliable"] for r in rows]
    for key, lab in (("R_fit_mm", "circle fit"), ("R_avg_mm", "arc length / bend")):
        axs[0].plot(p, [r[key] for r in rows], "o-", label=lab)
    for i, good in enumerate(ok):
        if not good:
            axs[0].plot(p[i], rows[i]["R_fit_mm"], "x", color="k")
    for label, _, size in REFERENCE_OBJECTS:
        axs[0].axhline(size / 2, ls="--", lw=0.8, color="0.5")
        axs[0].annotate(label, (p[-1], size / 2), fontsize=7, color="0.4", ha="right", va="bottom")
    axs[0].set(xlabel="chamber pressure (kPa)", ylabel="turn radius (mm)", yscale="log",
               title="turn radius (x: fit unreliable, < 90° bend)")
    axs[0].legend(fontsize=8)
    axs[1].plot(p, [r["circularity_rms_pct"] for r in rows], "o-")
    axs[1].set(xlabel="chamber pressure (kPa)", ylabel="rms off-circle (% of R)",
               title="circularity")
    for r in rows:
        prof = r["kappa_profile"]
        axs[2].plot([q["x_mm"] for q in prof], [1000 * q["kappa_per_mm"] for q in prof], "o-",
                    ms=3, label=f"{r['pressure_kPa']:g} kPa")
    axs[2].set(xlabel="undeformed position along finger (mm)", ylabel="curvature (1/m)",
               title="curvature along the finger (30 mm windows)")
    axs[2].legend(fontsize=7)
    for ax in axs:
        ax.grid(alpha=0.3)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    summary = out_dir / "turn_radius.png"
    fig.savefig(summary, dpi=140)
    plt.close(fig)
    return [sides, summary]
