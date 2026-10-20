"""Off-screen PyVista renders for checking a run by eye.

All images go to <run_dir>/renders/:
- mesh_tags.png       parts, cavity, clamp, tip, and bond interface, plus a
                      mid-plane cut-away so the chambers are visible
- deformed_XXXpct.png deformed shape at 25/50/75/100% of max pressure, side view
                      in the bending plane, fixed camera and color range, with
                      the undeformed outline in gray
- section_max.png     mid-plane section of the deformed state at max pressure,
                      colored by max principal stretch
- sweep_montage.png   the deformed frames in one image, plus sweep.gif
- contact_gaps_max.png slot walls between chambers at max pressure, colored by
                      the signed gap to the facing wall (red = interpenetration)

On request (render CLI --center-sections), renders/center/ holds sections on
the exact centre plane at every saved pressure level, with GIFs (and MP4s when
imageio-ffmpeg is installed); see render_center_sections.

Half models (mesh symmetry "half_y") are mirrored to the full finger first,
except for the centre sections, which use the half mesh's symmetry face.

Rendering uses the corner nodes only (linear tets), which is enough to see the
shape. Side views look along +y at the half y > mid-plane, so x runs left to
right from the clamp to the tip and the cut face is toward the camera.
"""

from pathlib import Path

import numpy as np

from ..mesh.fe_mesh import FEMesh, mirror_half, principal_stretches, tet10_jacobians
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
    """Distinct loaded levels nearest 25/50/75/100% of max, in order.

    A sweep with fewer than four loaded levels gets fewer frames rather than a
    repeated one.
    """
    p = np.asarray(pressures, float)
    loaded = np.flatnonzero(p > 0)
    idx = [int(loaded[np.argmin(np.abs(p[loaded] - f * p.max()))]) for f in FRACTIONS]
    return sorted(set(idx))


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


# --- centre-plane sections ----------------------------------------------------

CENTER_SIZES = {"full": (1600, 896), "zoom": (1600, 800)}   # multiples of 16 for MP4
CONTACT_COLORS = {"contact": "#00e5ff", "near": "#7fff00"}


def _text_box(actor):
    """White backing behind a text overlay so it stays readable over the mesh."""
    prop = actor.GetTextProperty()
    prop.SetBackgroundColor(1.0, 1.0, 1.0)
    prop.SetBackgroundOpacity(0.85)


def _rigid_xz(ref, cur):
    """Best-fit rotation R (2x2) and shift t with cur ~ R ref + t in the (x, z) plane."""
    P, Q = ref[:, [0, 2]], cur[:, [0, 2]]
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    if np.linalg.det(Vt.T @ U.T) < 0:
        Vt[-1] *= -1
    R = Vt.T @ U.T
    return R, qc - R @ pc


def _plane_camera(pl, center_xz, y, scale, up_xz=(0.0, 1.0)):
    """Parallel camera looking +y at (x, z) = center_xz, with an in-plane up direction."""
    c = (center_xz[0], y, center_xz[1])
    pl.camera_position = [(c[0], y - 500, c[2]), c, (up_xz[0], 0.0, up_xz[1])]
    pl.camera.parallel_projection = True
    pl.camera.parallel_scale = scale


def _facet_owners(tris, tets, n_nodes):
    """Index of the tet that each facet is a face of (matched on corners)."""
    n = np.int64(n_nodes)
    f = np.sort(tets[:, :4][:, [[0, 1, 3], [1, 2, 3], [2, 0, 3], [0, 2, 1]]], axis=2)
    f = f.reshape(-1, 3).astype(np.int64)
    keys = (f[:, 0] * n + f[:, 1]) * n + f[:, 2]
    t = np.sort(tris[:, :3], axis=1).astype(np.int64)
    q = (t[:, 0] * n + t[:, 1]) * n + t[:, 2]
    order = np.argsort(keys)
    hit = order[np.minimum(np.searchsorted(keys, q, sorter=order), len(order) - 1)]
    if not np.array_equal(keys[hit], q):
        raise ValueError("symmetry facets do not all match a tet face")
    return hit // 4


def _plane_section(mesh, tris, ids, points, lam, label):
    """On-plane tri6 facets as linear sub-triangles at `points`, plus element edges.

    `ids` are the facets' nodes; the polydata uses them as its points, in that
    order. Each sub-triangle carries the stretch of the tet that owns its facet.
    """
    pv = _pv()
    local = np.full(len(mesh.nodes), -1)
    local[ids] = np.arange(len(ids))
    owner = _facet_owners(tris, mesh.all_tets(), len(mesh.nodes))
    sub = local[np.vstack([tris[:, [0, 3, 5]], tris[:, [3, 1, 4]],
                           tris[:, [5, 4, 2]], tris[:, [3, 4, 5]]])]
    poly = pv.PolyData(points[ids], np.hstack([np.full((len(sub), 1), 3), sub]).ravel())
    poly.cell_data[label] = np.tile(lam[owner], 4)
    # tri6 edges as corner - mid-node - corner, each segment once
    seg = np.sort(local[tris[:, [0, 3, 3, 1, 1, 4, 4, 2, 2, 5, 5, 0]].reshape(-1, 2)], axis=1)
    seg = np.unique(seg, axis=0)
    edges = pv.PolyData(points[ids], lines=np.hstack([np.full((len(seg), 1), 2), seg]).ravel())
    return poly, edges


def _contour_section(mesh, y0, points, lam, label):
    """Full-model fallback: the y_ref = y0 isosurface of the deformed volume."""
    grid = volume_grid(mesh, points)
    grid.cell_data[label] = lam
    return grid.cell_data_to_point_data().contour([y0], scalars="y_ref"), None


def _zoom_x_range(mesh, zoom_slots, pad_mm=1.0):
    """Reference x range of the named slot walls, or the middle third of the finger."""
    names = [f"{s}_{side}" for s in zoom_slots for side in ("lo", "hi")
             if len(mesh.surfaces.get(f"{s}_{side}", []))]
    if names:
        x = mesh.nodes[np.unique(np.concatenate([mesh.surfaces[n].ravel() for n in names])), 0]
        return x.min() - pad_mm, x.max() + pad_mm, [s for s in zoom_slots
                                                     if f"{s}_lo" in names]
    x0, w = mesh.nodes[:, 0].min(), np.ptp(mesh.nodes[:, 0])
    return x0 + w / 3, x0 + 2 * w / 3, []


def _gap_summary(mesh, gap, on_plane, slots, tol):
    """Per slot: min signed gap over its wall nodes on the plane (None if none there)."""
    out = {}
    for s in slots:
        ids = np.intersect1d(np.unique(np.concatenate(
            [mesh.surfaces[f"{s}_lo"].ravel(), mesh.surfaces[f"{s}_hi"].ravel()])), on_plane)
        g = gap[ids]
        g = g[np.isfinite(g)]
        out[s] = (float(g.min()), int((g < tol).sum())) if len(g) else None
    return out


def render_center_sections(run_dir, zoom_slots=("slot04", "slot05", "slot06", "slot07"),
                           clim=None, gap_tol_mm=0.02, near_mm=0.1, frame_s=0.8, hold_s=2.5):
    """Sections on the exact centre plane at every saved pressure level, as frames and GIFs.

    Half models use the symmetry face of the half mesh itself (surface
    "symmetry", on y = plane_y_mm), so there is no contouring and no mirror
    seam. Full models fall back to the y_ref = mid-plane isosurface. Frames are
    the solved levels only; nothing is interpolated. Writes to
    <run_dir>/renders/center/:

    - center_full_KK_<p>kPa.png  whole finger, one fixed camera for all levels
    - center_zoom_KK_<p>kPa.png  `zoom_slots` with element edges, tracking the
                                 slots at a fixed zoom
    - center_full.gif, center_zoom.gif (and .mp4 with imageio-ffmpeg)

    Colored by max principal stretch on one color range for all frames (`clim`,
    default 1.0 to the largest value rounded up to 0.1). The undeformed section
    outline is gray. Slot-wall nodes on the plane are marked where the signed
    gap to the facing wall (contact.wall_gap_field) is below `gap_tol_mm`
    (contact) or `near_mm` (near).
    """
    pv = _pv()
    run_dir = Path(run_dir)
    out = run_dir / "renders" / "center"
    out.mkdir(parents=True, exist_ok=True)
    mesh = FEMesh.load(run_dir / "mesh" / "mesh.npz")
    d = np.load(run_dir / "displacements.npz")
    p, u = d["pressure_kPa"], d["u"]
    sym = mesh.info.get("symmetry")
    label = "max principal stretch"
    tets = mesh.all_tets()
    if sym:
        y0 = float(sym["plane_y_mm"])
        tris = mesh.surfaces["symmetry"]
        on_plane = np.unique(tris)
        how = "half-model symmetry face"
    else:
        y0 = _mid_y(mesh)
        on_plane = np.flatnonzero(np.abs(mesh.nodes[:, 1] - y0) < 1e-6)
        how = "isosurface of the full model"
    slots = [s for s in mesh.info.get("slots", []) if len(mesh.surfaces.get(f"{s}_lo", []))]
    x_lo, x_hi, zoom_found = _zoom_x_range(mesh, zoom_slots)

    frames = []
    for k in range(len(p)):
        pts = mesh.nodes + u[k]
        lam = principal_stretches(mesh.nodes, pts, tets)[:, -1]
        sec, edges = (_plane_section(mesh, tris, on_plane, pts, lam, label) if sym
                      else _contour_section(mesh, y0, pts, lam, label))
        gap = wall_gap_field(mesh, pts)
        g = gap[on_plane]
        frames.append({
            "p": float(p[k]), "pts": pts, "sec": sec, "edges": edges,
            "contact": on_plane[g < gap_tol_mm],
            "near": on_plane[(g >= gap_tol_mm) & (g < near_mm)],
            "gaps": _gap_summary(mesh, gap, on_plane, slots, gap_tol_mm),
            "inverted": int((tet10_jacobians(pts, tets).min(1) <= 0).sum()),
        })
    if clim is None:
        top = max(float(np.max(f["sec"][label])) for f in frames)
        clim = (1.0, max(1.1, np.ceil(top * 10) / 10))

    ref_poly = (_plane_section(mesh, tris, on_plane, mesh.nodes, np.ones(len(tets)), label)[0]
                if sym else _contour_section(mesh, y0, mesh.nodes, np.ones(len(tets)), label)[0])
    outline = ref_poly.extract_feature_edges(boundary_edges=True, feature_edges=False,
                                             manifold_edges=False, non_manifold_edges=False)
    xmin, xmax, zmin, zmax = _xz_bounds(np.vstack([ref_poly.points]
                                                  + [f["sec"].points for f in frames]))
    full_bounds = (xmin, xmax, zmin, zmax + 0.08 * (zmax - zmin))   # headroom for the text
    # zoom: the on-plane wall nodes of the zoom slots (nodes near the plane in
    # their x range for full models). Each level gets their best-fit rigid motion
    # from the reference; the zoom camera turns with it and the gray outline is
    # moved by it, at one zoom for all levels.
    walls = [np.intersect1d(np.unique(np.concatenate(
        [mesh.surfaces[f"{s}_{side}"].ravel() for side in ("lo", "hi")])), on_plane)
        for s in zoom_found]
    zoom_ids = (np.concatenate(walls) if walls else np.flatnonzero(
        (mesh.nodes[:, 0] >= x_lo) & (mesh.nodes[:, 0] <= x_hi)
        & (np.abs(mesh.nodes[:, 1] - y0) < 1.0)))
    aspect = CENTER_SIZES["zoom"][0] / CENTER_SIZES["zoom"][1]
    fits, need = [], []
    for f in frames:
        R, t = _rigid_xz(mesh.nodes[zoom_ids], f["pts"][zoom_ids])
        local = (f["pts"][zoom_ids][:, [0, 2]] - t) @ R      # R^T (x - t): slot frame
        lo, hi = local.min(0), local.max(0)
        fits.append((R, t, (lo + hi) / 2))
        w, h = 1.1 * (hi - lo) + 2.0
        need.append(max(h, w / (0.82 * aspect)))
    zoom_scale = 0.56 * max(need)

    levels = ", ".join(f"{q:g}" for q in p)
    zoom_name = (f"{zoom_found[0]}-{zoom_found[-1]}" if zoom_found
                 else f"x_ref {x_lo:.0f}-{x_hi:.0f} mm")
    written = {"full": [], "zoom": []}
    for k, f in enumerate(frames):
        contact_txt = (f"slot-wall nodes on the plane: {len(f['contact'])} in contact "
                       f"(gap < {gap_tol_mm:g} mm), {len(f['near'])} near (< {near_mm:g} mm)")
        head = (f"p = {f['p']:g} kPa, exact centre plane y = {y0:g} mm ({how})\n"
                f"solved level {k + 1} of {len(p)} ({levels} kPa), not interpolated\n")
        if f["inverted"]:
            head += f"WARNING: {f['inverted']} inverted elements in the model\n"
        R, t, c_local = fits[k]
        # box centre at 41% of the width, clear of the color bar
        c_local = c_local + [0.09 * 2 * zoom_scale * aspect, 0.0]
        turn = round(float(np.degrees(np.arctan2(R[1, 0], R[0, 0]))))
        zoom_gaps = "; ".join(
            f"{s}: " + (f"min {f['gaps'][s][0]:.3f} mm" if f["gaps"].get(s) else "none")
            for s in zoom_found)
        zoom_outline = outline.copy()
        moved = np.array(outline.points)
        moved[:, [0, 2]] = moved[:, [0, 2]] @ R.T + t
        zoom_outline.points = moved
        views = {
            "full": (outline, head + contact_txt, False),
            "zoom": (zoom_outline,
                     head + f"zoom on {zoom_name}, element edges shown, view turned "
                     f"{turn} deg with the slots\n" + contact_txt
                     + (f"\non-plane wall gap {zoom_gaps}" if zoom_gaps else ""),
                     True),
        }
        for view, (ref_line, text, show_edges) in views.items():
            size = CENTER_SIZES[view]
            pl = pv.Plotter(off_screen=True, window_size=size)
            # at y = y0, so push the outline behind and edges and markers in front
            pl.add_mesh(ref_line.translate((0, 0.5, 0)), color="#9e9e9e", line_width=1.5,
                        opacity=0.6)
            pl.add_mesh(f["sec"], scalars=label, clim=clim, cmap="magma",
                        show_edges=show_edges and f["edges"] is None, edge_color="#999999",
                        scalar_bar_args={"title": label, **BAR})
            if show_edges and f["edges"] is not None:
                pl.add_mesh(f["edges"].translate((0, -0.05, 0)), color="#999999",
                            line_width=0.8)
            for kind in ("near", "contact"):
                if len(f[kind]):
                    # plain points: sphere points do not render off-screen here
                    pl.add_points(f["pts"][f[kind]] + [0, -0.2, 0],
                                  color=CONTACT_COLORS[kind],
                                  point_size=10 if view == "zoom" else 6)
            if view == "zoom":
                for s, ids in zip(zoom_found, walls):
                    top = ids[mesh.nodes[ids, 2] >= mesh.nodes[ids, 2].max() - 1e-6]
                    pl.add_point_labels(f["pts"][top].mean(0, keepdims=True) + [0, -0.3, 0],
                                        [s], font_size=14, show_points=False,
                                        shape_opacity=0.6, always_visible=True)
            _text_box(pl.add_text(text, font_size=10 if view == "full" else 9))
            _text_box(pl.add_text(
                "gray: undeformed outline" + (", rigidly moved with the slots" if show_edges
                                              else "")
                + "   cyan: on-plane slot-wall nodes in contact   green: near",
                position="lower_left", font_size=8))
            if view == "full":
                _side(pl, full_bounds, aspect=size[0] / size[1])
            else:
                _plane_camera(pl, R @ c_local + t, y0, zoom_scale, up_xz=R[:, 1])
            path = out / f"center_{view}_{k:02d}_{f['p']:g}kPa.png"
            pl.screenshot(str(path))
            pl.close()
            written[view].append(path)

    paths = written["full"] + written["zoom"]
    for view, pngs in written.items():
        paths.append(_gif(pngs, out / f"center_{view}.gif", frame_s, hold_s))
        mp4 = _mp4(pngs, out / f"center_{view}.mp4", frame_s, hold_s)
        if mp4:
            paths.append(mp4)
    summary = [{"pressure_kPa": f["p"], "contact_nodes": len(f["contact"]),
                "near_nodes": len(f["near"]), "inverted": f["inverted"],
                "slot_min_gap_on_plane_mm": {s: (v[0] if v else None)
                                             for s, v in f["gaps"].items()}}
               for f in frames]
    return paths, summary


def _montage(frames, path):
    import matplotlib.image as mpimg
    imgs = [mpimg.imread(str(f))[..., :3] for f in frames]
    while len(imgs) % 2:
        imgs.append(np.ones_like(imgs[0]))
    rows = [np.hstack(imgs[i:i + 2]) for i in range(0, len(imgs), 2)]
    mpimg.imsave(str(path), np.vstack(rows))


def _gif(frames, path, frame_s=0.8, hold_s=None):
    """Looping GIF with Pillow, frame_s per frame and hold_s on the last one."""
    from PIL import Image
    imgs = [Image.open(f).convert("RGB") for f in frames]
    ms = [round(1000 * frame_s)] * len(imgs)
    ms[-1] = round(1000 * (hold_s or frame_s))
    imgs[0].save(str(path), save_all=True, append_images=imgs[1:], duration=ms, loop=0)
    return path


def _mp4(frames, path, frame_s=0.8, hold_s=None, fps=10):
    """MP4 at `fps` with each frame repeated for its duration; None without imageio-ffmpeg."""
    try:
        import imageio.v2 as imageio
        import imageio_ffmpeg  # noqa: F401
    except ImportError:
        return None
    reps = [round(frame_s * fps)] * len(frames)
    reps[-1] = round((hold_s or frame_s) * fps)
    with imageio.get_writer(str(path), fps=fps, macro_block_size=16) as w:
        for f, n in zip(frames, reps):
            img = imageio.imread(str(f))[..., :3]
            for _ in range(n):
                w.append_data(img)
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
