"""Design metrics from nodal displacements, independent of the solver.

Definitions (finger axis `a`, e.g. +x; bend side `b`, e.g. -z, the base side):

- Tip position: area-weighted centroid of the distal tip face, deformed.
- Tip normal: area-weighted outward normal of the deformed tip face
  (initially +a). It is the end tangent of the finger's neutral line.
- Bend angle: rotation of the tip normal in the bending plane spanned by
  a and b, theta = atan2(n.b, n.a), in degrees. 0 is straight, 90 means the
  tip face points toward the base side, and values past 180 mean the finger
  has curled back on itself (angles are unwrapped across pressure steps).
  Positive is toward the base, which is the expected direction.
- Out-of-plane angle: asin of the tip normal's component on a x b, a
  check for twist or sideways bending.
- Cavity volume: divergence theorem over the deformed cavity facets.
- Peak wall stretch: max principal stretch at tet10 centroids, per part.
"""

import numpy as np

from ..mesh.fe_mesh import cavity_volume, principal_stretches, tri6_to_linear

AXIS_VECTORS = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}


def axis_vector(spec):
    """'+x', '-z', 'x' -> unit vector."""
    sign = -1.0 if spec.startswith("-") else 1.0
    return sign * np.array(AXIS_VECTORS[spec.lstrip("+-")], float)


def face_frame(points, tris6):
    """Area-weighted centroid and unit normal of a triangulated face."""
    t = points[tri6_to_linear(tris6)]
    n = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]) / 2
    area = np.linalg.norm(n, axis=1)
    centroid = (t.mean(1) * area[:, None]).sum(0) / area.sum()
    normal = n.sum(0)
    return centroid, normal / np.linalg.norm(normal)


def bend_angle_deg(normal, axis, bend):
    return float(np.degrees(np.arctan2(normal @ bend, normal @ axis)))


def step_metrics(mesh, u, axis="x", bend_toward="-z"):
    """Metrics for one displacement field u (N, 3) mm."""
    a, b = axis_vector(axis), axis_vector(bend_toward)
    x = mesh.nodes + u
    c0, _ = face_frame(mesh.nodes, mesh.surfaces["tip"])
    c, n = face_frame(x, mesh.surfaces["tip"])
    vol, _ = cavity_volume(x, mesh.surfaces["cavity"])
    stretch = {}
    for name, tets in mesh.parts.items():
        lam = principal_stretches(mesh.nodes, x, tets)
        stretch[name] = {"max": float(lam[:, -1].max()), "min": float(lam[:, 0].min())}
    return {
        "tip_position_mm": c.tolist(),
        "tip_displacement_mm": (c - c0).tolist(),
        "tip_normal": n.tolist(),
        "bend_angle_deg": bend_angle_deg(n, a, b),
        "out_of_plane_deg": float(np.degrees(np.arcsin(np.clip(n @ np.cross(a, b), -1, 1)))),
        "cavity_volume_mm3": float(vol),
        "principal_stretch": stretch,
    }


def unwrap_angles(steps):
    """Unwrap bend angles in place across a sweep."""
    ang = np.unwrap(np.radians([s["bend_angle_deg"] for s in steps]))
    for s, v in zip(steps, np.degrees(ang)):
        s["bend_angle_deg"] = float(v)
    return steps


def synthetic_bend(nodes, angle_deg, axis_min_mm=0.0, length_mm=None, neutral_z_mm=0.0):
    """Constant-curvature bend toward -z about the y axis, for testing.

    The section at x = axis_min stays fixed and the end at axis_min + length
    rotates by `angle_deg`. The neutral surface is z = neutral_z.
    Returns displacement (N, 3).
    """
    x, z = nodes[:, 0] - axis_min_mm, nodes[:, 2] - neutral_z_mm
    L = length_mm or float(x.max())
    theta_end = np.radians(angle_deg)
    if abs(theta_end) < 1e-12:
        return np.zeros_like(nodes)
    R = L / theta_end
    th = x / R
    r = R + z
    new = nodes.copy()
    new[:, 0] = axis_min_mm + r * np.sin(th)
    new[:, 2] = neutral_z_mm - R + r * np.cos(th)
    return new - nodes
