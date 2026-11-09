"""Split an existing cut holding vector by its own mean-stress nodal force.

The complement includes body/load/support conventions, not just deviatoric
stress. These are diagnostic force paths, with no new mesh or material law.
"""
import numpy as np
from verify_high_order import prepare
from spatial_release_control import fractions


def partition(mesh, active, stress, holding, order):
    active = np.asarray(active, dtype=bool)
    stress = np.asarray(stress, dtype=float)
    holding = np.asarray(holding, dtype=float)
    if (active.shape != (len(mesh['tetrahedra']),) or
            stress.shape != (27*int(active.sum()), 6) or
            holding.shape != (3*len(mesh['nodes']),) or
            not np.isfinite(stress).all() or not np.isfinite(holding).all()):
        raise ValueError('invalid existing full field or holding vector')
    if order not in ('mean-first', 'complement-first'):
        raise ValueError('unknown stress-force release order')
    quadrature, _ = prepare(mesh)
    mean_force = np.zeros_like(holding)
    position = 0
    for element in np.flatnonzero(active):
        q = quadrature[element]
        isotropic = np.zeros((27, 6))
        isotropic[:, :3] = stress[position:position+27, :3].mean(axis=1)[:, None]
        local = -np.einsum('qai,qa,q->i', q['b'], isotropic, q['weight'])
        dofs = (q['nodes'][:, None]*3+np.arange(3)).ravel()
        np.add.at(mean_force, dofs, local)
        position += 27
    complement = holding-mean_force
    return ((mean_force, complement) if order == 'mean-first'
            else (complement, mean_force))


def load(final, holding, first, second, s):
    a, b = fractions(s)
    final = np.asarray(final)
    if s == 0.:
        return final+np.asarray(holding)
    if s == 1.:
        return final.copy()
    return final+(1-a)*np.asarray(first)+(1-b)*np.asarray(second)
