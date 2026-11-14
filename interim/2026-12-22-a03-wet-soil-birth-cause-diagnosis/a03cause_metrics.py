"""Fixed C-order observations. These metrics never authorize acceptance."""
import math
import numpy as np


def flat(value):
    a = np.asarray(value, dtype=float)
    if not np.isfinite(a).all():
        raise ValueError('nonfinite field')
    return a.ravel(order='C')


def dot(a, b):
    a, b = flat(a), flat(b)
    if a.shape != b.shape:
        raise ValueError('component correspondence differs')
    return math.fsum(float(x)*float(y) for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def metric(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError('field correspondence differs')
    d = flat(a)-flat(b)
    raw = norm(b)
    denominator = max(raw, 1e-12)
    absolute = norm(d)
    return dict(absolute_l2=absolute, reference_l2=raw, denominator=denominator,
                relative_l2=absolute/denominator, maximum_absolute=float(max(abs(d), default=0.)),
                components=int(d.size), denominator_floor_active=raw < 1e-12)


def compare(left, right, categories):
    if len(left) != len(right) or not left:
        raise ValueError('missing accepted state')
    keys = sorted(left[0][1])
    for (ra, fa), (rb, fb) in zip(left, right):
        if (ra['stage'], ra['progress']) != (rb['stage'], rb['progress']):
            raise ValueError('accepted coordinates differ; interpolation is forbidden')
        if sorted(fa) != keys or sorted(fb) != keys:
            raise ValueError('missing field')
    result = {}
    for key in keys:
        a = np.stack([f[key] for _, f in left]); b = np.stack([f[key] for _, f in right])
        if a.shape != b.shape:
            raise ValueError('field correspondence differs')
        if key in categories:
            result[key] = dict(mismatches=int(np.count_nonzero(a != b)), components=int(a.size))
        else:
            result[key] = metric(a, b)
    return dict(states=len(left), fields=result, numerical_acceptance=False)


def common(rows, steps):
    result = [item for item in rows if item[0]['progress']*steps == int(item[0]['progress']*steps)]
    if not result:
        raise ValueError('no common coordinates')
    return result
