"""Observation-only fixed-order arithmetic; no production or legacy imports."""
import math
from pathlib import Path
import numpy as np


def dot(a, b):
    a = np.asarray(a, dtype=float).ravel(order='C')
    b = np.asarray(b, dtype=float).ravel(order='C')
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('finite, identically shaped vectors required')
    return math.fsum(float(x)*float(y) for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def mv(a, b):
    return np.array([dot(row, b) for row in np.asarray(a)])


def cross(a):
    x, y, z = a
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def frame(points):
    p = np.asarray(points)
    n = np.cross(p[1]-p[0], p[2]-p[0])
    return n/norm(n)


def arm_observation(points, e, t):
    n = frame(points)
    return dict(normal=n.tolist(), e_m=e, thickness_m=t, arm_m=(e*n).tolist(),
                finite_within_original_range=bool(math.isfinite(e) and abs(e) <= t/2),
                original_shared_cutoff_m=64*np.finfo(float).eps*t)


def shared_guard(arms, thicknesses):
    """Source-level observation of one common ID in declaration order."""
    previous = None
    comparisons = []
    for arm, t in zip(arms, thicknesses):
        arm = np.asarray(arm)
        if previous is not None:
            distance = norm(previous-arm)
            cutoff = 64*np.finfo(float).eps*t
            comparisons.append(dict(distance_m=distance, cutoff_m=cutoff,
                                    refused=bool(distance > cutoff)))
            if distance > cutoff:
                return dict(accepted=False, comparisons=comparisons)
        previous = arm
    return dict(accepted=True, comparisons=comparisons)


def rigid_basis(positions):
    p = np.asarray(positions)
    result = np.zeros((6*len(p), 6))
    for i, x in enumerate(p):
        result[6*i:6*i+3, :3] = np.eye(3)
        result[6*i:6*i+3, 3:] = -cross(x)
        result[6*i+3:6*i+6, 3:] = np.eye(3)
    return result


def wrench(positions, vector):
    f = np.asarray(vector).reshape(-1, 6)
    forces = [math.fsum(float(x) for x in f[:, a]) for a in range(3)]
    moments = [math.fsum(float(np.cross(x, row[:3])[a]+row[3+a])
                         for x, row in zip(positions, f)) for a in range(3)]
    return dict(force_kn=forces, moment_kn_m=moments)


def document_text(path):
    return Path(path).read_text(encoding='utf-8').replace('\r\n', '\n')


def observed_order(values):
    if len(values) < 3:
        return None
    a, b, c = values[-3:]
    d1, d2 = abs(b-a), abs(c-b)
    return math.log2(d1/d2) if d1 > 0 and d2 > 0 else None
