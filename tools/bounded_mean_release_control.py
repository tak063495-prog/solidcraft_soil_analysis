"""Mean-force-guided partition that can only withdraw original holding loads.

Clip the force partition, never material states. Both vectors retain each
original component's sign, and sum to the original holding vector. Original
zero components stay zero at every intermediate coordinate.
"""
import numpy as np
from mean_stress_release_control import partition as raw_partition, fractions, load


def partition(mesh, active, stress, holding, order):
    if order not in ('mean-first', 'complement-first'):
        raise ValueError('unknown bounded stress-force release order')
    mean, _ = raw_partition(mesh, active, stress, holding, 'mean-first')
    h = np.asarray(holding)
    mean = np.minimum(np.maximum(mean, np.minimum(h, 0.)), np.maximum(h, 0.))
    complement = h-mean
    return (mean, complement) if order == 'mean-first' else (complement, mean)
