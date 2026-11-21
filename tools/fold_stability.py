"""Read-only diagnostics for a nonsymmetric equilibrium path R(u,lambda)=0.

SVD always uses K itself. Coordinate scales and the moving tangent plane are
reported, so an augmented condition number is never compared across hidden
changes of units. These diagnostics do not change a constitutive law or a
solver tolerance, and a local fold says nothing about disconnected equilibria.
"""
import numpy as np

FORCE_TOL = 1e-8
LINEAR_TOL = 1e-10


def array(value, shape=None):
    a = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(a)) or (shape is not None and a.shape != shape):
        raise ValueError('nonfinite or incorrectly shaped fold data')
    return a


def spectrum(a):
    s = np.linalg.svd(a, compute_uv=False)
    return dict(singular_values=s.tolist(), minimum=float(s[-1]), maximum=float(s[0]),
                condition_number=float(s[0]/s[-1]) if s[-1] else None,
                numerical_rank=int(np.sum(s > max(a.shape)*np.finfo(float).eps*s[0])))


def point(k, h, previous_tangent=None, displacement_scale_m=1., load_scale=1.,
          force_scale_kn=1., contributions=()):
    """Physical left/right weak modes and scaled bordered-system diagnostics.

    z=(u/displacement_scale_m, lambda/load_scale), ||dz/ds||=1.
    previous_tangent is dz/ds, with its sign fixed by the preceding path point.
    The bottom row is the actual moving plane used for correction; the null
    direction is independently reconstructed from [K,H], without inverting K.
    """
    k = array(k)
    if k.ndim != 2 or k.shape[0] != k.shape[1] or not 1 <= len(k) <= 256:
        raise ValueError('fold SVD is limited to square matrices of 1..256 DOF')
    n = len(k); h = array(h, (n,))
    scales = array([displacement_scale_m, load_scale, force_scale_kn])
    if np.any(scales <= 0) or np.linalg.norm(h) == 0:
        raise ValueError('positive scales and a nonzero load direction are required')
    left, singular, right = np.linalg.svd(k, full_matrices=True)
    v = right[-1].copy(); w = left[:, -1].copy()
    # Deterministic paired sign: K v = sigma w.
    if v[np.argmax(np.abs(v))] < 0: v = -v; w = -w
    force_block = np.column_stack((k*displacement_scale_m, h*load_scale))/force_scale_kn
    _, _, vh = np.linalg.svd(force_block, full_matrices=True)
    tangent = vh[-1]
    if previous_tangent is None:
        if abs(tangent[-1]) < 64*np.finfo(float).eps:
            raise ValueError('orientation at a fold needs a previous branch direction')
        if tangent[-1] < 0: tangent = -tangent
        plane = tangent.copy()
    else:
        plane = array(previous_tangent, (n+1,))
        if abs(np.linalg.norm(plane)-1) > 1e-10:
            raise ValueError('moving-plane direction must have unit scaled norm')
        alignment = float(tangent@plane)
        if abs(alignment) < 64*np.finfo(float).eps:
            raise ValueError('ambiguous branch orientation')
        if alignment < 0: tangent = -tangent
    augmented = np.vstack((force_block, plane))
    pieces = []
    total = np.zeros_like(k)
    for item in contributions:
        local = array(item['matrix'], k.shape); total += local
        pieces.append({key: val for key, val in item.items() if key != 'matrix'} | dict(
            w_k_v=float(w@local@v), k_v_norm=float(np.linalg.norm(local@v)),
            v_k_v=float(v@local@v)))
    return dict(dimension=n, tangent_kind='actual_nonsymmetric',
                stiffness=spectrum(k), right_weak_mode=v.tolist(), left_weak_mode=w.tolist(),
                right_mode_residual=float(np.linalg.norm(k@v-singular[-1]*w)),
                left_mode_residual=float(np.linalg.norm(k.T@w-singular[-1]*v)),
                w_dot_h=float(w@h), normalized_transversality=float(abs(w@h)/np.linalg.norm(h)),
                stiffness_asymmetry=float(np.linalg.norm(k-k.T)/max(np.linalg.norm(k), 1e-300)),
                scaled_tangent=tangent.tolist(), dlambda_ds=float(load_scale*tangent[-1]),
                tangent_force_block_residual=float(np.linalg.norm(force_block@tangent)),
                moving_plane=plane.tolist(), augmented=spectrum(augmented),
                displacement_scale_m=float(displacement_scale_m), load_scale=float(load_scale),
                force_scale_kn=float(force_scale_kn), contributions=pieces,
                contribution_sum_relative_error=(float(np.linalg.norm(total-k)/max(np.linalg.norm(k),1e-300))
                                                 if pieces else None))


def audit_load(actual, final, holding, lam):
    """Check every physical component, including restrained/inactive DOFs."""
    actual = array(actual); final = array(final, actual.shape); holding = array(holding, actual.shape)
    if actual.ndim != 1 or not np.isfinite(lam): raise ValueError('invalid affine load')
    expected = final+(1-lam)*holding
    error = actual-expected
    scale = max(float(np.linalg.norm(expected)), 1.)
    if np.max(np.abs(error), initial=0.) > FORCE_TOL*scale or np.linalg.norm(error) > FORCE_TOL*scale:
        raise ValueError('component load violates Ffinal+(1-lambda)H')
    return dict(maximum_component_error_kn=float(np.max(np.abs(error), initial=0.)),
                relative_full_error=float(np.linalg.norm(error)/scale),
                remaining_fraction=float(1-lam))


def judge(points):
    """Evidence gate. Missing paths and failed solves are never fold evidence."""
    if len(points) < 3:
        return dict(conclusion='undetermined', reason='fewer than three accepted points', brackets=[])
    for p in points:
        for key, limit in [('relative_force_residual', FORCE_TOL),
                           ('force_block_linear_residual', LINEAR_TOL),
                           ('augmented_linear_residual', LINEAR_TOL)]:
            value = p.get(key)
            if value is None or not np.isfinite(value) or value > limit or value < 0:
                return dict(conclusion='undetermined', reason='missing or failed '+key, brackets=[])
        if not p.get('accepted', False):
            return dict(conclusion='undetermined', reason='rejected point in accepted path', brackets=[])
    brackets = []
    for i, (a,b) in enumerate(zip(points, points[1:])):
        if a['diagnostic']['dlambda_ds'] > 0 > b['diagnostic']['dlambda_ds']:
            brackets.append(dict(indices=[i,i+1], lambdas=[a['lambda'], b['lambda']]))
    if not brackets:
        return dict(conclusion='undetermined', reason='no oriented dlambda/ds sign change', brackets=[])
    for bracket in brackets:
        lo,hi=bracket['indices']
        for p in points[max(0,lo-1):min(len(points),hi+2)]:
            d=p['diagnostic']; n=d['dimension']
            eps=64*np.finfo(float).eps
            if (d['augmented']['numerical_rank'] != n+1 or
                    d['normalized_transversality'] <= eps or
                    d['stiffness']['numerical_rank'] < n-1):
                return dict(conclusion='undetermined',reason='degenerate or nontransverse bracket',brackets=brackets)
    return dict(conclusion='local_fold_bracket', brackets=brackets,
                reason='accepted oriented tangent changes + to - with transverse load and regular border',
                limitation='requires refinement/independent agreement; no global nonexistence or dynamic stability claim')


def distribution(a):
    a=array(a)
    if not a.size: return dict(count=0,minimum=None,maximum=None,quantiles=None)
    return dict(count=int(a.size), minimum=float(a.min()), maximum=float(a.max()),
                quantiles=dict(zip(['0','0.01','0.1','0.5','0.9','0.99','1'],
                                   np.quantile(a,[0,.01,.1,.5,.9,.99,1]).tolist())))


def low_stress(points):
    """All GP values, including frozen removed points; thresholds are diagnostic."""
    records=[]
    for p in points:
        stress=array(p['stress_kpa'], (6,)); mean=float(np.mean(stress[:3]))
        e=float(p['void_ratio']); kappa=float(p['kappa'])
        if not np.isfinite(e) or e <= -1 or not np.isfinite(kappa) or kappa <= 0:
            raise ValueError('invalid MCC bulk-modulus parameters')
        records.append(dict(p, mean_effective_pressure_kpa=mean,
                            bulk_modulus_kpa=(1+e)*mean/kappa))
    active=[p for p in records if p['active']]
    result=dict(schema='soil-fem-low-stress-diagnostic/1',diagnostic_only=True,
                formula='K=(1+e)*p_prime/kappa; no pressure floor',points=records)
    for name, values in [('all',records),('active',active)]:
        result[name]=dict(pressure_kpa=distribution([p['mean_effective_pressure_kpa'] for p in values]),
                         bulk_modulus_kpa=distribution([p['bulk_modulus_kpa'] for p in values]),
                         near_zero_pressure=[dict(threshold_kpa=t,count=sum(p['mean_effective_pressure_kpa']<=t for p in values),
                                                  point_ids=[p['point'] for p in values if p['mean_effective_pressure_kpa']<=t])
                                             for t in (1.,.1,.01,.001,1e-6)])
    candidates=[p for p in active if p['bulk_modulus_kpa']>0 and 'weak_mode_strain_norm' in p]
    if len(candidates)>1:
        x=np.log([p['bulk_modulus_kpa'] for p in candidates]); y=np.array([p['weak_mode_strain_norm'] for p in candidates])
        correlation=float(np.corrcoef(x,y)[0,1]) if np.std(x)>0 and np.std(y)>0 else None
    else:correlation=None
    result['pearson_log_bulk_vs_weak_mode_strain']=correlation
    result['interpretation']='correlation is diagnostic, not proof that K tends to zero or of stability'
    return result
