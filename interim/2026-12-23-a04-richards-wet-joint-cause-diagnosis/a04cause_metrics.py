"""Observation only: fixed state / ID / component order, no acceptance override."""
import json
import math
from pathlib import Path
import numpy as np
CATEGORIES = {'cell_ids', 'face_keys', 'active_elements', 'material_return', 'joint_status', 'joint_generation', 'joint_active', 'body_active', 'ties', 'contact_status', 'contact_open', 'joint_tied', 'generation', 'archived_generation_count', 'anchor_0_taut', 'frozen_status', 'exterior_faces', 'wet_exterior_faces'}


def flat(a):
    a = np.asarray(a, dtype=float)
    if not np.isfinite(a).all():
        raise ValueError('nonfinite observation')
    return a.ravel(order='C')


def norm(a):
    return math.sqrt(math.fsum(float(x)*float(x) for x in flat(a)))


def metric(a, b):
    if np.shape(a) != np.shape(b):
        raise ValueError('component correspondence differs')
    a, b = flat(a), flat(b)
    delta = a-b
    raw = norm(b)
    absolute = norm(delta)
    return dict(absolute_l2=absolute, reference_l2=raw, denominator=max(raw, 1e-12),
                relative_l2=absolute/max(raw, 1e-12), maximum_absolute=float(max(abs(delta), default=0.)),
                components=len(a), denominator_floor_active=raw < 1e-12)


def rows(path):
    with Path(path).open(encoding='utf-8') as f:
        return [(r['row'], r['fields']) for r in map(json.loads, f)]


def coordinate(row):
    return json.dumps(row, sort_keys=True, separators=(',', ':'), allow_nan=False)


def common(series, grid):
    wanted = [coordinate(r) for r, _ in grid]
    lookup = {coordinate(r): (r, f) for r, f in series}
    if len(lookup) != len(series) or any(k not in lookup for k in wanted):
        raise ValueError('duplicate or missing physical coordinate; no interpolation')
    return [lookup[k] for k in wanted]


def compare(a, b, keys=None, categories=CATEGORIES):
    if not a or len(a) != len(b):
        raise ValueError('missing accepted states')
    keys = sorted(a[0][1]) if keys is None else sorted(keys)
    for (ra, fa), (rb, fb) in zip(a, b):
        if coordinate(ra) != coordinate(rb):
            raise ValueError('accepted coordinates differ')
        if keys == sorted(a[0][1]) and (sorted(fa) != keys or sorted(fb) != keys):
            raise ValueError('field inventory differs')
        for key in keys:
            if key not in fa or key not in fb or np.shape(fa[key]) != np.shape(fb[key]):
                raise ValueError('field/ID correspondence differs: '+key)
    result = {}
    for key in keys:
        left = np.concatenate([np.asarray(f[key]).ravel(order='C') for _, f in a])
        right = np.concatenate([np.asarray(f[key]).ravel(order='C') for _, f in b])
        result[key] = dict(mismatches=int(np.count_nonzero(left != right)), components=len(left)) if key in categories else metric(left, right)
    return dict(states=len(a), fields=result, numerical_acceptance=False)


def convergence(series, grid, keys):
    selected = [common(s, grid) for s in series]
    pairs = [compare(selected[i], selected[i+1], keys) for i in range(2)]
    return dict(common_coordinates=[r for r, _ in grid], interpolation=False, extrapolation=False,
                pairs=pairs, orders={k: math.log2(pairs[0]['fields'][k]['absolute_l2']/pairs[1]['fields'][k]['absolute_l2'])
                if pairs[0]['fields'][k]['absolute_l2'] > 0 and pairs[1]['fields'][k]['absolute_l2'] > 0 else None for k in keys},
                numerical_acceptance=False)
