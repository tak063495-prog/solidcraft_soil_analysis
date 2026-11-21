"""Numerical ordering of an existing holding vector; no new physical load."""
import numpy as np


def fractions(s):
    if not np.isfinite(s) or not 0 <= s <= 1:
        raise ValueError('release coordinate outside [0,1]')
    return min(2*s, 1.), max(2*s-1., 0.)


def partition(nodes, holding, strategy):
    xyz = np.asarray(nodes, dtype=float)
    h = np.asarray(holding, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or h.shape != (xyz.size,) or not np.isfinite(xyz).all() or not np.isfinite(h).all():
        raise ValueError('invalid nodal holding vector')
    if strategy == 'left-first':
        chosen = xyz[:, 0] < (xyz[:, 0].min()+xyz[:, 0].max())/2
    elif strategy == 'center-first':
        xy = xyz[:, :2]
        center = (xy.min(axis=0)+xy.max(axis=0))/2
        radius = np.max(np.ptp(xy, axis=0))*.51
        chosen = np.linalg.norm(xy-center, axis=1) < radius
    else:
        raise ValueError('unknown spatial release strategy')
    mask = np.repeat(chosen, 3)
    first = np.where(mask, h, 0.)
    second = np.where(mask, 0., h)
    if not np.any(first) or not np.any(second):
        raise ValueError('both holding groups must be nonzero')
    return first, second, np.flatnonzero(chosen).tolist()


def load(final, first, second, s):
    a, b = fractions(s)
    return np.asarray(final)+(1-a)*np.asarray(first)+(1-b)*np.asarray(second)
