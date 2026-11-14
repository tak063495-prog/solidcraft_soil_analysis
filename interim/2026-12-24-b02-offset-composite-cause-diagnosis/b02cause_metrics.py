"""Fixed coordinate/ID/component order and fsum; never relax an old gate."""
import math
from pathlib import Path
import numpy as np


def dot(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('finite arrays with identical shapes required')
    return math.fsum(float(x)*float(y) for x, y in zip(a.ravel(order='C'), b.ravel(order='C'), strict=True))


def norm(a):
    return math.sqrt(dot(a, a))


def metric(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    dot(a, b)
    delta = a-b
    raw = norm(b)
    return dict(absolute_l2=norm(delta), maximum_absolute=float(np.max(np.abs(delta), initial=0.)),
                reference_l2=raw, denominator=max(raw, 1e-12), denominator_floor_active=raw < 1e-12,
                relative_l2=norm(delta)/max(raw, 1e-12), components=a.size)


def common(rows, steps=4):
    result = {}
    for row in rows:
        key = (row['stage_index'], row['progress'])
        if row['progress'] not in tuple(j/steps for j in range(1, steps+1)):
            continue
        if key in result:
            raise ValueError('duplicate physical path coordinate')
        result[key] = row
    if len(result) != 5*steps:
        raise ValueError('missing common path coordinates; no interpolation')
    return result


def document_text(path):
    # Text comparison contracts use universal newlines, including CRLF checkout.
    return Path(path).read_text(encoding='utf-8-sig')
