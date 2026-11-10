"""Stress-path diagnostics; no stability verdict or full-removal acceptance.

Compression-positive MCC invariants and signs follow the existing associated
flow and logarithmic pc hardening law. Native fields are observations only.
"""
import numpy as np


def unpack(states, kind):
    if kind == 'native':
        stress = np.asarray([s['effective_stress_compression_positive_kpa'] for s in states], dtype=float)
        pc = np.asarray([s['preconsolidation_pressure_kpa'] for s in states], dtype=float)
        void = np.asarray([s['void_ratio'] for s in states], dtype=float)
        plastic = np.asarray([s['plastic_strain_compression_positive'] for s in states], dtype=float)
    elif kind == 'independent':
        a = np.asarray(states, dtype=float).reshape(-1, 14)
        stress, pc, void, plastic = a[:, :6], np.exp(a[:, 6]), np.expm1(a[:, 7]), a[:, 8:14]
    else:
        raise ValueError('unknown material map')
    if (stress.shape != (len(pc), 6) or plastic.shape != stress.shape
            or not all(np.isfinite(a).all() for a in (stress, pc, void, plastic))
            or np.any(pc <= 0) or np.any(void <= 0)):
        raise ValueError('nonfinite or inadmissible material observation')
    return stress, pc, void, plastic


def invariants(states, kind, critical_state_m):
    stress, pc, void, plastic = unpack(states, kind)
    m = np.asarray(critical_state_m, dtype=float)
    if not np.isfinite(m).all() or np.any(m <= 0):
        raise ValueError('invalid input MCC M')
    p = stress[:, :3].mean(axis=1)
    dev = stress.copy(); dev[:, :3] -= p[:, None]
    q2 = 1.5*(np.sum(dev[:, :3]**2, axis=1)+2*np.sum(dev[:, 3:]**2, axis=1))
    # For f=q²/M²+p(p-pc), tr(df/dsigma)=2p-pc.
    # d ln(pc)=(1+e)/(lambda-kappa) d eps_v^p.
    return dict(p=p, q=np.sqrt(q2), pc=pc, void=void, plastic=plastic,
                p_over_pc=p/pc, flow_volume_normal=2*p-pc,
                yield_over_pc_squared=(q2/m**2+p*(p-pc))/pc**2)


def metrics(a):
    a = np.asarray(a)
    if not a.size or not np.isfinite(a).all():
        raise ValueError('empty or nonfinite diagnostic field')
    return dict(minimum=float(a.min()), maximum=float(a.max()), l2=float(np.linalg.norm(a)))


def describe(states, kind, critical_state_m, active_gp, weights, previous=None):
    f = invariants(states, kind, critical_state_m)
    pi = np.asarray(active_gp, dtype=int)
    w = np.asarray(weights, dtype=float)
    if (len(pi) != len(w) or len(set(pi.tolist())) != len(pi)
            or np.any(pi < 0) or np.any(pi >= len(f['p']))
            or not np.isfinite(w).all() or np.any(w <= 0)):
        raise ValueError('invalid observed active GP inventory')
    dry = f['flow_volume_normal'][pi] < 0
    result = dict(active_points=len(pi), dry_side_points=int(dry.sum()),
                  dry_side_volume_fraction=float(w[dry].sum()/w.sum()),
                  invariants={key: metrics(f[key][pi]) for key in
                              ('p', 'q', 'pc', 'void', 'p_over_pc', 'flow_volume_normal', 'yield_over_pc_squared')},
                  physical_failure_confirmed=False, numerical_acceptance=False)
    if previous is not None:
        before = invariants(previous, kind, critical_state_m)
        depv = (f['plastic'][pi, :3]-before['plastic'][pi, :3]).sum(axis=1)
        dlpc = np.log(f['pc'][pi]/before['pc'][pi])
        result['increment'] = dict(plastic_volume_strain=metrics(depv), log_pc=metrics(dlpc),
            dilating_points=int(np.sum(depv < 0)), softening_points=int(np.sum(dlpc < 0)),
            dilating_and_softening_points=int(np.sum((depv < 0) & (dlpc < 0))),
            signed_dilation_pc_product=metrics(depv*dlpc))
    return result


def continuum_terms(states, kind, parameters):
    """Reevaluate existing MCC rate equations at observed states, not a solver.

    No branch selection or discrete-return-map tangent is implied. Both the
    elastic and hypothetical plastic-loading tensors are reported separately.
    """
    f = invariants(states, kind, parameters['critical_state_m'])
    stress, pc, void, _ = unpack(states, kind)
    lam, kap, nu, m = (np.asarray(parameters[k], dtype=float) for k in
                      ('lambda', 'kappa', 'poisson', 'critical_state_m'))
    lam, kap, nu, m = (np.broadcast_to(a, (len(pc),)) for a in (lam, kap, nu, m))
    if (not all(np.isfinite(a).all() for a in (lam, kap, nu, m))
            or np.any(kap <= 0) or np.any(lam <= kap) or np.any(nu <= -1) or np.any(nu >= .5)
            or np.any(f['p'] <= 0)):
        raise ValueError('inadmissible MCC rate diagnostic input')
    bulk = f['p']*(1+void)/kap
    shear = 3*bulk*(1-2*nu)/(2*(1+nu))
    elastic = np.zeros((len(pc), 6, 6))
    elastic[:, :3, :3] = (bulk-2*shear/3)[:, None, None]
    for i in range(3):
        elastic[:, i, i] += 2*shear
        elastic[:, i+3, i+3] = shear
    dev = stress.copy(); dev[:, :3] -= f['p'][:, None]
    normal = 3*dev/m[:, None]**2
    normal[:, :3] += f['flow_volume_normal'][:, None]/3
    normal[:, 3:] *= 2
    dn = np.einsum('qij,qj->qi', elastic, normal)
    ndn = np.einsum('qi,qi->q', normal, dn)
    hardening = f['p']*pc*(1+void)/(lam-kap)*f['flow_volume_normal']
    denominator = ndn+hardening
    if np.any(denominator == 0):
        raise ValueError('zero consistency denominator; no tangent invented')
    plastic = elastic-np.einsum('qi,qj,q->qij', dn, dn, 1/denominator)
    return dict(elastic=elastic, plastic_loading=plastic, normal=normal,
                hardening=hardening, denominator=denominator,
                normal_elastic_normal=ndn)
