"""Explicit, provenance-selected correction; never reinterpret old saved runs."""
import numpy as np
import fem_reference as fem
from hydrostatic_cut_volume_reference import force


def correction(mesh, active, quadrature, level, gamma, soil=None):
    result = np.zeros((len(mesh['nodes']), 3))
    if level is None:
        return result
    xyz = np.asarray(mesh['nodes'])
    for ei in np.flatnonzero(active):
        if soil is not None and not soil[ei]:
            continue
        q = quadrature[ei]; points = xyz[q['nodes']]; vertices = points[:4]
        if not vertices[:, 2].min() < level < vertices[:, 2].max():
            continue
        if len(points) == 10:
            mids = np.asarray([(vertices[i]+vertices[j])/2 for i, j in fem.EDGES])
            scale = np.linalg.norm(vertices[1:]-vertices[0], axis=1).max()
            if np.any(np.linalg.norm(points[4:]-mids, axis=1) > 32*np.finfo(float).eps*scale):
                continue
        p = gamma*np.maximum(level-q['position'][:, 2], 0)
        sampled = -np.einsum('qi,q,q->i', q['b'][:, :3].sum(axis=1), p, q['weight'])
        result[q['nodes']] += force(mesh, ei, level, gamma)['force']-sampled.reshape(-1, 3)
    return result


def enabled(run):
    value = run.get('hydrostatic_volume_integration')
    assert value in (None, 'affine_cut_partial'), 'unknown hydrostatic volume integration'
    return value is not None
