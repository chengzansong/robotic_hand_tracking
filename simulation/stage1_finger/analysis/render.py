"""Off-screen PyVista renders for checking a run by eye.

All images go to <run_dir>/renders/:
- mesh_tags.png       parts, cavity, clamp, tip, and bond interface, plus a
                      mid-plane cut-away so the chambers are visible
- deformed_XXXpct.png deformed shape at 25/50/75/100% of max pressure, side view
                      in the bending plane, fixed camera and color range, with
                      the undeformed outline in gray
- section_max.png     mid-plane section of the deformed state at max pressure,
                      colored by max principal stretch
- sweep_montage.png   the deformed frames in one image (plus sweep.gif when
                      imageio is installed)
- contact_gaps_max.png slot walls between chambers at max pressure, colored by
                      the signed gap to the facing wall (red = interpenetration)

Half models (mesh symmetry "half_y") are mirrored to the full finger first.

Rendering uses the corner nodes only (linear tets), which is enough to see the
shape. Side views look along +y at the half y > mid-plane, so x runs left to
right from the clamp to the tip and the cut face is toward the camera.
"""

from pathlib import Path

import numpy as np

from ..mesh.fe_mesh import FEMesh, mirror_half, principal_stretches
from .contact import wall_gap_field

PART_COLORS = ["#9ecae1", "#bdbdbd", "#c7e9c0", "#fdd0a2"]
TAG_COLORS = {"cavity": "#f16913", "clamp": "#cb181d", "tip": "#238b45",
              "interface": "#ffd92f", "slot": "#6a51a3"}
FRACTIONS = (0.25, 0.5, 0.75, 1.0)
EDGES = {"show_edges": True, "edge_color": "#555555", "line_width": 0.3}
BAR = {"vertical": True, "position_x": 0.88, "position_y": 0.15, "height": 0.7}


def _pv():
    import pyvista as pv
    pv.OFF_SCREEN = True
    return pv


def volume_grid(mesh, points=None):
    """Linear-tet UnstructuredGrid with part id and reference y as data."""
    pv = _pv()
    tets = mesh.all_tets()[:, :4]
    cells = np.hstack([np.full((len(tets), 1), 4), tets]).ravel()
    pts = mesh.nodes if points is None else points
    grid = pv.UnstructuredGrid(cells, np.full(len(tets), pv.CellType.TETRA), pts)
    grid.cell_data["part"] = np.concatenate(
        [np.full(len(v), i) for i, v in enumerate(mesh.parts.values())])
    grid.point_data["y_ref"] = mesh.nodes[:, 1]
    return grid


def surface_poly(mesh, name, points=None):
    pv = _pv()
    tris = mesh.surfaces[name][:, :3]
    faces = np.hstack([np.full((len(tris), 1), 3), tris]).ravel()
    poly = pv.PolyData(mesh.nodes if points is None else points, faces)
    poly.point_data["y_ref"] = mesh.nodes[:, 1]
    return poly


def _skin(grid):
    return grid.extract_surface(algorithm="dataset_surface")


def _outline(grid):
    """Sharp feature edges only, for a faint undeformed silhouette."""
    return _skin(grid).extract_feature_edges(
        feature_angle=45, boundary_edges=False, non_manifold_edges=False,
        feature_edges=True, manifold_edges=False)


def _mid_y(mesh):
    return float(mesh.nodes[:, 1].min() + np.ptp(mesh.nodes[:, 1]) / 2)


def _far_half(data, y_mid):
    """Keep reference y > y_mid (works on deformed data too)."""
    return data.clip_scalar(scalars="y_ref", value=y_mid, invert=False)


def _add_parts(pl, mesh, grid, legend=False, **kw):
    for i, name in enumerate(mesh.parts):
        sub = grid.extract_cells(np.flatnonzero(grid.cell_data["part"] == i))
        if sub.n_cells:
            pl.add_mesh(sub, color=PART_COLORS[i % len(PART_COLORS)],
                        label=name if legend else None, **kw)


def _iso(pl, mesh, from_tip=True):
    c = mesh.nodes.mean(0)
    span = np.ptp(mesh.nodes, axis=0).max()
    sx = 1 if from_tip else -1
    pl.camera_position = [tuple(c + span * np.array([0.9 * sx, -0.9, 0.7])), tuple(c), (0, 0, 1)]


def _side(pl, bounds, aspect):
    """Fixed parallel camera looking +y, covering bounds (xmin, xmax, zmin, zmax)."""
    xmin, xmax, zmin, zmax = bounds
    c = ((xmin + xmax) / 2, 7.5, (zmin + zmax) / 2)
    pl.camera_position = [(c[0], c[1] - 500, c[2]), c, (0, 0, 1)]
    pl.camera.parallel_projection = True
    pl.camera.parallel_scale = 0.56 * max(zmax - zmin, (xmax - xmin) / (0.82 * aspect))


def _xz_bounds(points):
    return (points[:, 0].min(), points[:, 0].max(), points[:, 2].min(), points[:, 2].max())


def render_mesh_tags(mesh, path):
    """Tagged mesh from both ends, cavity alone, and mid-plane cut-aways."""
    pv = _pv()
    grid = volume_grid(mesh)
    y_mid = _mid_y(mesh)
    cut = _far_half(grid, y_mid)
    cav_cut = _far_half(surface_poly(mesh, "cavity"), y_mid)
    pl = pv.Plotter(off_screen=True, shape=(2, 3), window_size=(2400, 1200), border=False)

    for col, from_tip, tag in ((0, True, "tip"), (1, False, "clamp")):
        pl.subplot(0, col)
        pl.add_text(f"From the {'tip' if from_tip else 'clamp'} end: top (blue), base (gray), "
                    f"{tag} ({'green' if tag == 'tip' else 'red'})", font_size=10)
        _add_parts(pl, mesh, _skin(grid), legend=col == 0, **EDGES)
        # nudge the tag face outward along the finger axis so it does not z-fight
        shift = np.zeros(3)
        shift["xyz".index(mesh.info.get("fem_tags", {}).get("axis", "x"))] = 0.05
        pl.add_mesh(surface_poly(mesh, tag).translate(shift if from_tip else -shift),
                    color=TAG_COLORS[tag])
        if col == 0:
            pl.add_legend(bcolor="white", size=(0.2, 0.12), loc="upper right")
        _iso(pl, mesh, from_tip)

    pl.subplot(0, 2)
    pl.add_text("Cavity (orange), bond (yellow), contact slot walls (purple)", font_size=10)
    pl.add_mesh(_skin(grid), color="#dddddd", opacity=0.12)
    pl.add_mesh(surface_poly(mesh, "cavity"), color=TAG_COLORS["cavity"])
    if len(mesh.surfaces.get("interface", [])):
        pl.add_mesh(surface_poly(mesh, "interface"), color=TAG_COLORS["interface"], opacity=0.6)
    for name in _slot_surfaces(mesh):
        pl.add_mesh(surface_poly(mesh, name), color=TAG_COLORS["slot"])
    _iso(pl, mesh)

    pl.subplot(1, 0)
    pl.add_text(f"Cut at y = {y_mid:g} mm, cavity walls orange", font_size=10)
    _add_parts(pl, mesh, cut, **EDGES)
    pl.add_mesh(cav_cut, color=TAG_COLORS["cavity"])
    _iso(pl, mesh)

    pl.subplot(1, 1)
    pl.add_text("Cut, side view: clamp left, tip right", font_size=10)
    _add_parts(pl, mesh, cut, **EDGES)
    pl.add_mesh(cav_cut, color=TAG_COLORS["cavity"])
    _side(pl, _xz_bounds(mesh.nodes), aspect=1.0)

    pl.subplot(1, 2)
    pl.add_text("Cut, side view, no edges", font_size=10)
    _add_parts(pl, mesh, cut)
    pl.add_mesh(cav_cut, color=TAG_COLORS["cavity"])
    pl.add_mesh(_far_half(surface_poly(mesh, "clamp"), y_mid), color=TAG_COLORS["clamp"])
    pl.add_mesh(_far_half(surface_poly(mesh, "tip"), y_mid), color=TAG_COLORS["tip"])
    _side(pl, _xz_bounds(mesh.nodes), aspect=1.0)

    pl.screenshot(str(path))
    pl.close()
    return path


def _slot_surfaces(mesh):
    return [f"{s}_{side}" for s in mesh.info.get("slots", []) for side in ("lo", "hi")
            if len(mesh.surfaces.get(f"{s}_{side}", []))]


def render_contact(mesh, pressure_kPa, u_max, path, clim=(-0.2, 1.0)):
    """Contact footprint at max pressure: one gap map per slot.

    Each panel is the -x wall of one slot drawn in its reference (y, z)
    coordinates and colored by the signed gap to the facing wall in the
    deformed state. Blue is open, near zero (yellow) is touching, red is
    interpenetration. The walls are edge-on in the side views, hence a map.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    slots = [s for s in mesh.info.get("slots", []) if len(mesh.surfaces.get(f"{s}_lo", []))]
    if not slots:
        return None
    gap = wall_gap_field(mesh, mesh.nodes + u_max)
    cols = 5
    rows = -(-len(slots) // cols)
    fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 3.0 * rows + 0.6), squeeze=False)
    im = None
    for ax, slot in zip(axes.flat, slots):
        tris = mesh.surfaces[f"{slot}_lo"][:, :3]
        ids = np.unique(tris)
        local = np.full(len(mesh.nodes), -1)
        local[ids] = np.arange(len(ids))
        y, z = mesh.nodes[ids, 1], mesh.nodes[ids, 2]
        g = gap[ids]
        tri = mtri.Triangulation(y, z, local[tris])
        im = ax.tripcolor(tri, np.nan_to_num(g, nan=clim[1]), cmap="RdYlBu",
                          vmin=clim[0], vmax=clim[1], shading="gouraud")
        ax.set_title(f"{slot}: min {np.nanmin(g):.3f} mm", fontsize=9)
        ax.set_aspect("equal")
        ax.set_xlabel("y (mm)", fontsize=8)
        ax.set_ylabel("z (mm)", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in list(axes.flat)[len(slots):]:
        ax.axis("off")
    fig.suptitle(f"Slot-wall gap at p = {pressure_kPa:.3g} kPa (the -x wall of each slot, "
                 f"reference coordinates). Red = interpenetration.", fontsize=11)
    fig.colorbar(im, ax=axes, shrink=0.8, label="signed gap (mm)")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _frame_indices(pressures):
    p = np.asarray(pressures)
    return [int(np.argmin(np.abs(p - f * p.max()))) for f in FRACTIONS]


def render_sweep(mesh, pressures_kPa, u, renders_dir, color_by="displacement"):
    """Deformed side views at 25/50/75/100% of max pressure plus a montage.

    `u` is (steps, N, 3) nodal displacement in mm, matching `pressures_kPa`.
    """
    pv = _pv()
    renders_dir = Path(renders_dir)
    outline = _outline(volume_grid(mesh))
    idx = _frame_indices(pressures_kPa)
    p_max = float(np.max(pressures_kPa))
    bounds = _xz_bounds(np.vstack([mesh.nodes] + [mesh.nodes + u[i] for i in idx]))
    size = (1400, 800)

    tets = mesh.all_tets()
    if color_by == "stretch":
        scal = [principal_stretches(mesh.nodes, mesh.nodes + u[i], tets)[:, -1] for i in idx]
        label = "max principal stretch"
    else:
        scal = [np.linalg.norm(u[i], axis=1) for i in idx]
        label = "|u| (mm)"
    clim = (min(s.min() for s in scal), max(s.max() for s in scal))

    frames = []
    for i, s in zip(idx, scal):
        grid = volume_grid(mesh, mesh.nodes + u[i])
        (grid.cell_data if color_by == "stretch" else grid.point_data)[label] = s
        pl = pv.Plotter(off_screen=True, window_size=size)
        pl.add_mesh(outline, color="#9e9e9e", line_width=1.5, opacity=0.6)
        pl.add_mesh(_skin(grid), scalars=label, clim=clim, cmap="viridis",
                    scalar_bar_args={"title": label, **BAR})
        pct = 100 * pressures_kPa[i] / p_max
        pl.add_text(f"p = {pressures_kPa[i]:.3g} kPa ({pct:.0f}% of max)", font_size=12)
        _side(pl, bounds, aspect=size[0] / size[1])
        out = renders_dir / f"deformed_{round(pct):03d}pct.png"
        pl.screenshot(str(out))
        pl.close()
        frames.append(out)
    _montage(frames, renders_dir / "sweep_montage.png")
    _gif(frames, renders_dir / "sweep.gif")
    return frames


def render_section(mesh, pressure_kPa, u_max, path):
    """Mid-plane section at max pressure: chamber ballooning and wall thinning."""
    pv = _pv()
    y_mid = _mid_y(mesh)
    label = "max principal stretch"
    # A true planar section in reference coordinates: the y_ref = y_mid isosurface.
    ref = volume_grid(mesh).contour([y_mid], scalars="y_ref")
    grid = volume_grid(mesh, mesh.nodes + u_max)
    grid.cell_data[label] = principal_stretches(mesh.nodes, mesh.nodes + u_max,
                                                mesh.all_tets())[:, -1]
    section = grid.cell_data_to_point_data().contour([y_mid], scalars="y_ref")
    size = (1600, 900)
    pl = pv.Plotter(off_screen=True, window_size=size)
    pl.add_mesh(ref.extract_feature_edges(boundary_edges=True, feature_edges=False,
                                          manifold_edges=False, non_manifold_edges=False),
                color="#9e9e9e", line_width=1.5, opacity=0.7)
    pl.add_mesh(section, scalars=label, cmap="magma",
                scalar_bar_args={"title": label, **BAR})
    pl.add_text(f"Mid-plane section y = {y_mid:g} mm at p = {pressure_kPa:.3g} kPa",
                font_size=12)
    _side(pl, _xz_bounds(np.vstack([mesh.nodes, mesh.nodes + u_max])), size[0] / size[1])
    pl.screenshot(str(path))
    pl.close()
    return path


def _montage(frames, path):
    import matplotlib.image as mpimg
    imgs = [mpimg.imread(str(f))[..., :3] for f in frames]
    while len(imgs) % 2:
        imgs.append(np.ones_like(imgs[0]))
    rows = [np.hstack(imgs[i:i + 2]) for i in range(0, len(imgs), 2)]
    mpimg.imsave(str(path), np.vstack(rows))


def _gif(frames, path):
    try:
        import imageio.v2 as imageio
    except ImportError:
        return None
    imageio.mimsave(str(path), [imageio.imread(str(f)) for f in frames], duration=0.8)
    return path


def render_run(run_dir):
    """(Re)generate every render for a run folder. Returns the written paths."""
    run_dir = Path(run_dir)
    renders = run_dir / "renders"
    renders.mkdir(exist_ok=True)
    half = FEMesh.load(run_dir / "mesh" / "mesh.npz")
    mesh = mirror_half(half)
    written = [render_mesh_tags(mesh, renders / "mesh_tags.png")]
    disp = run_dir / "displacements.npz"
    if disp.exists():
        d = np.load(disp)
        p, u = d["pressure_kPa"], d["u"]
        if half.info.get("symmetry"):
            _, u = mirror_half(half, u)
        if len(p) and p.max() > 0:
            written += render_sweep(mesh, p, u, renders)
            k = int(np.argmax(p))
            written.append(render_section(mesh, p[k], u[k], renders / "section_max.png"))
            written.append(renders / "sweep_montage.png")
            c = render_contact(mesh, p[k], u[k], renders / "contact_gaps_max.png")
            if c:
                written.append(c)
    return written
