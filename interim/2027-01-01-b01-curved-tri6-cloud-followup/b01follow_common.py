"""Fixed-order observations, separate from LAPACK solve/eigenvalue arithmetic."""
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'c6b1176b1147b8f242ef979f6c4a6a70e6829579'
EDGES = ((0, 1), (1, 2), (2, 0))
SECTION = dict(id='b01follow', young_kpa=3e7, poisson=.2, thickness_m=.1,
               unit_weight_kn_m3=0., drilling_factor=1e-4)


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def dot(a, b):
    a, b = np.asarray(a, float).ravel(order='C'), np.asarray(b, float).ravel(order='C')
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('finite equally sized operands required')
    return math.fsum(float(x)*float(y) for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def mm(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[0]:
        raise ValueError('matrix dimensions differ')
    return np.array([[dot(row, col) for col in b.T] for row in a])


def mv(a, b):
    return np.array([dot(row, b) for row in np.asarray(a)])


def metric(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape:
        raise ValueError('shape mismatch')
    absolute, reference = norm(a-b), norm(b)
    return dict(absolute_l2=absolute, reference_l2=reference,
                relative_l2=absolute/max(reference, 1e-12),
                maximum_absolute=float(np.max(np.abs(a-b), initial=0.)))


def rotate(v, q):
    value = np.asarray(v, float)
    return mm(value.reshape(-1, 3), np.asarray(q).T).reshape(value.shape)


def rigid(x):
    x = np.asarray(x); center = [math.fsum(x[:, a])/len(x) for a in range(3)]
    modes = []
    for axis in np.eye(3):
        modes.append(np.tile(np.r_[axis, [0., 0., 0.]], (len(x), 1)).ravel())
    for axis in np.eye(3):
        modes.append(np.column_stack((np.cross(np.tile(axis, (len(x), 1)), x-center),
                                      np.tile(axis, (len(x), 1)))).ravel())
    return np.column_stack(modes)
