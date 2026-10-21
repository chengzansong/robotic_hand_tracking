import numpy as np

def finger_curvature(points):
    """
    Estimate constant curvature from 4 MediaPipe 3D landmarks.

    points: NumPy array with shape (4, 3), in meters, ordered MCP -> tip.
    Returns: unsigned curvature in 1/meters.
    """
    points = np.asarray(points, dtype=float)

    # center points to (0,0,0) since their absolute position doesnt matter
    centroid = np.mean(points, axis=0)
    centered = points - centroid

    # find best fit plane in 3d
    _, singular_values, vh = np.linalg.svd(centered)

    if singular_values[0] < 1e-8:
        return 0.0

    # project 3d points on plane
    basis = vh[:2].T
    xy = centered @ basis

    # arc length along the finger, chord from base to tip
    arc_length = float(np.sum(np.linalg.norm(np.diff(xy, axis=0), axis=1)))
    if arc_length < 1e-6:
        return 0.0

    chord = float(np.linalg.norm(xy[-1] - xy[0]))
    ratio = chord / arc_length

    # straight (or numerically longer chord than polyline)
    if ratio >= 1.0:
        return 0.0

    # c = L * sinc(theta/2)  =>  chord/arc = sin(x)/x with x = theta/2
    if ratio <= 0.0:
        half_theta = np.pi
    else:
        lo, hi = 0.0, np.pi
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if np.sin(mid) / mid > ratio:
                lo = mid
            else:
                hi = mid
        half_theta = 0.5 * (lo + hi)

    theta = 2.0 * half_theta
    return float(theta / arc_length)