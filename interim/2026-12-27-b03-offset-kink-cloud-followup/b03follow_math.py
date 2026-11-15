"""Observation-only deterministic scalar reductions. No production hooks."""
import json
import math
from pathlib import Path
import numpy as np


def serial(value):
    if isinstance(value, np.ndarray):
        return serial(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    return value


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='\n') as stream:
        stream.write(json.dumps(serial(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def dot(a, b):
    a = np.asarray(a, dtype=float).ravel(order='C')
    b = np.asarray(b, dtype=float).ravel(order='C')
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('finite equal-shape vectors required')
    return math.fsum(float(x)*float(y) for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def mv(a, b):
    return np.array([dot(row, b) for row in np.asarray(a)])


def cross(a):
    x, y, z = a
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def arm_map(arm, nodes=6):
    a = np.eye(6*nodes)
    for i in range(nodes):
        a[6*i:6*i+3, 6*i+3:6*i+6] = cross(arm)
    return a


def rigid(positions):
    r = np.zeros((6*len(positions), 6))
    for i, x in enumerate(positions):
        r[6*i:6*i+3, :3] = np.eye(3)
        r[6*i:6*i+3, 3:] = -cross(x)
        r[6*i+3:6*i+6, 3:] = np.eye(3)
    return r


def assemble(parts, size):
    terms = {}
    for dofs, k in parts:
        for i, gi in enumerate(dofs):
            for j, gj in enumerate(dofs):
                terms.setdefault((gi, gj), []).append(float(k[i, j]))
    result = np.zeros((size, size))
    for (i, j), values in sorted(terms.items()):
        result[i, j] = math.fsum(values)
    return result


def scatter(parts, size):
    terms = [[] for _ in range(size)]
    for dofs, values in parts:
        for i, value in zip(dofs, values, strict=True):
            terms[i].append(float(value))
    return np.array([math.fsum(values) for values in terms])


def guard(arms, thicknesses):
    previous = None
    comparisons = []
    for arm, thickness in zip(arms, thicknesses, strict=True):
        arm = np.asarray(arm)
        if previous is not None:
            delta = norm(arm-previous)
            cutoff = 64*np.finfo(float).eps*thickness
            refused = delta > cutoff
            comparisons.append(dict(distance_m=delta, cutoff_m=cutoff, refused=refused))
            if refused:
                return dict(accepted=False, comparisons=comparisons)
        previous = arm
    return dict(accepted=True, comparisons=comparisons)


def document_text(path):
    return Path(path).read_text(encoding='utf-8').replace('\r\n', '\n')
