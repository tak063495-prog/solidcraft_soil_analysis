"""Diagnostic LP column generation for GP equilibrium, never a continuum FELA certificate.

Only numpy/scipy already locked by this repository are used. Frozen-pc MCC is
an ellipsoid. Convex combinations of its support points give feasible stresses;
LP dual multipliers and analytic support functions bound the objective below.
"""
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix, hstack, vstack

P = np.array([1., 1., 1., 0., 0., 0.])/3
D = np.eye(6)
D[:3, :3] -= 1/3
QMAP = np.diag([1., 1., 1., np.sqrt(2), np.sqrt(2), np.sqrt(2)]) @ D * np.sqrt(1.5)


def ellipse_metric(m):
    if not np.isfinite(m) or m <= 0:
        raise ValueError('positive finite M required')
    return QMAP.T @ QMAP / m**2 + np.outer(P, P)


def ellipse_support(direction, pc, m):
    """Argmax and support of q²/M²+(p-pc/2)² <= (pc/2)²."""
    direction = np.asarray(direction, dtype=float)
    if direction.shape != (6,) or not np.all(np.isfinite(direction)) or not np.isfinite(pc) or pc <= 0:
        raise ValueError('finite six-vector and positive pc required')
    center = np.r_[np.full(3, pc/2), np.zeros(3)]
    v = np.linalg.solve(ellipse_metric(m), direction)
    length = np.sqrt(max(float(direction @ v), 0.))
    stress = center if length == 0 else center + pc/2 * v/length
    return stress, float(stress @ direction)


def cone_ray(direction, m):
    """Maximise direction:stress on the p=1 cross-section of q<=M*p."""
    direction = np.asarray(direction, dtype=float)
    dev = direction.copy(); dev[:3] -= direction[:3].mean()
    # q metric restricted to deviatoric stresses, with engineering shear duality.
    v = dev / np.array([1.5, 1.5, 1.5, 3., 3., 3.])
    length = np.sqrt(max(float(dev @ v), 0.))
    stress = np.r_[np.ones(3), np.zeros(3)]
    if length:
        stress += m*v/length
    return stress, float(stress @ direction)



def intersection_support(direction, pc, m):
    """Support of the frozen ellipse intersected with q<=M*p (closed apex).

    Below p=pc/2 the cone is active; above it the ellipse is active.
    With negative hydrostatic objective the optimum is the critical ring
    or the apex. With nonnegative objective it is the ellipse support.
    """
    if np.sum(direction[:3]) >= 0:
        return ellipse_support(direction, pc, m)
    ray, value = cone_ray(direction, m)
    stress = ray*(pc/2) if value > 0 else np.zeros(6)
    return stress, float(stress@direction)


def yield_values(stress, pc, m, kind):
    stress = np.asarray(stress); p = stress @ P; q = np.linalg.norm(stress @ QMAP.T, axis=-1)
    ellipse = (q/(m*pc))**2 + (p/pc)*(p/pc-1)
    if kind == 'ellipse':
        return ellipse
    if kind == 'intersection':
        return np.maximum(ellipse, (q-m*p)/(m*pc))
    return q-m*p


def discrete_bound(b, weights, final, holding, pc, m, kind='ellipse', max_iterations=160, gap_tolerance=2e-7):
    """min mu, -sum w B^T sigma=final+mu*holding, 0<=mu<=1.

    Returns a feasible restricted LP and an independently evaluated analytic
    dual bound. Finite numerical residuals are reported, never hidden in a
    certified physical bound. Each GP is independent; no traction continuity
    or pointwise divergence is enforced by this displacement-FE equilibrium.
    """
    b = np.asarray(b); weights = np.asarray(weights); final = np.asarray(final); holding = np.asarray(holding)
    pc = np.broadcast_to(pc, (len(b),)); m = np.broadcast_to(m, (len(b),))
    if kind not in ('ellipse', 'cone', 'intersection') or b.ndim != 3 or b.shape[1] != 6 or len(weights) != len(b):
        raise ValueError('invalid bound layout or kind')
    if not all(np.all(np.isfinite(x)) for x in (b, weights, final, holding, pc, m)) or np.any(weights <= 0) or np.any(pc <= 0) or np.any(m <= 0):
        raise ValueError('finite data and positive GP weights, pc, M required')
    nq, _, nd = b.shape
    if final.shape != (nd,) or holding.shape != (nd,):
        raise ValueError('load layout differs')
    a = -np.einsum('qai,q->qia', b, weights)
    support_function = intersection_support if kind == 'intersection' else ellipse_support
    columns, owners, stresses = [], [], []
    def add(q, stress):
        columns.append(a[q] @ stress); owners.append(q); stresses.append(stress)
    for q in range(nq):
        if kind != 'cone':
            add(q, np.r_[np.full(3, pc[q]/2), np.zeros(3)])
        else:
            add(q, np.r_[np.ones(3), np.zeros(3)])
        for d in np.vstack((np.eye(6), -np.eye(6))):
            stress, _ = support_function(d, pc[q], m[q]) if kind != 'cone' else cone_ray(d, m[q])
            add(q, stress)
    trace = []; result = None; lower = None
    for iteration in range(max_iterations):
        nc = len(columns)
        equilibrium = hstack([csc_matrix(np.asarray(columns).T), csc_matrix(-holding[:, None])], format='csc')
        if kind != 'cone':
            simplex = csc_matrix((np.ones(nc), (owners, np.arange(nc))), shape=(nq, nc))
            eq = vstack([equilibrium, hstack([simplex, csc_matrix((nq, 1))])], format='csc')
            rhs = np.r_[final, np.ones(nq)]
        else:
            eq = equilibrium; rhs = final
        result = linprog(np.r_[np.zeros(nc), 1.], A_eq=eq, b_eq=rhs,
                         bounds=[(0, None)]*nc+[(0, 1)], method='highs',
                         options={'dual_feasibility_tolerance': 1e-9, 'primal_feasibility_tolerance': 1e-9})
        if not result.success:
            # A restricted master can be infeasible even when the cone problem
            # is feasible; return that limitation instead of claiming collapse.
            return dict(status='restricted_lp_unresolved', solver_message=result.message,
                        iterations=iteration+1, trace=trace, physical_bound_certified=False)
        y = result.eqlin.marginals[:nd]
        support_sum = 0.; violations = []; additions = []
        for q in range(nq):
            d = a[q].T @ y
            stress, support = support_function(d, pc[q], m[q]) if kind != 'cone' else cone_ray(d, m[q])
            support_sum += support
            violation = support + (result.eqlin.marginals[nd+q] if kind != 'cone' else 0.)
            violations.append(violation)
            if violation > 1e-11:
                additions.append((q, stress))
        if kind != 'cone':
            lower = float(y@final-support_sum+min(0., 1.+y@holding))
        elif max(violations) <= 1e-11:
            lower = float(y@final+min(0., 1.+y@holding))
        trace.append(dict(iteration=iteration, columns=nc, mu_feasible=float(result.fun), mu_dual_lower=lower,
                          maximum_pricing_violation=float(max(violations))))
        if lower is not None and result.fun-lower <= gap_tolerance:
            break
        if not additions:
            break
        for q, stress in additions:
            add(q, stress)
    # result belongs to the columns BEFORE the final pricing additions.
    nc = len(result.x)-1
    stress = np.zeros((nq, 6))
    for q, value, s in zip(owners[:nc], result.x[:nc], stresses[:nc]):
        stress[q] += value*s
    mu = float(result.fun)
    residual = np.einsum('qia,qa->i', a, stress)-final-mu*holding
    return dict(status='bounded_discrete_model' if lower is not None and mu-lower <= gap_tolerance else 'feasible_discrete_only',
                kind=kind, mu_feasible=mu, mu_dual_lower=lower,
                lambda_LB=1-mu, lambda_optimum_upper=None if lower is None else 1-lower,
                objective_gap=None if lower is None else mu-lower,
                stress_kpa=stress.tolist(), dual_equilibrium=y.tolist(),
                equilibrium_max_abs_kn=float(np.max(np.abs(residual))),
                equilibrium_relative=float(np.linalg.norm(residual)/max(np.linalg.norm(final+mu*holding), 1.)),
                maximum_yield=float(np.max(yield_values(stress, pc, m, kind))),
                minimum_p_kpa=float(np.min(stress@P)), iterations=len(trace), trace=trace,
                physical_bound_certified=False,
                limitation='GP weak equilibrium only, not H(div) continuum equilibrium; frozen pc, not hardening MCC path')


def mechanism_bound(b, weights, mode, final, holding, pc, m):
    """Virtual work bound for the frozen-pc quadrature model, with its sign."""
    de = -np.einsum('qai,i->qa', np.asarray(b), np.asarray(mode))
    support = [ellipse_support(e, c, mm)[1] for e, c, mm in zip(de, pc, m)]
    dissipation = float(np.dot(weights, support))
    final_work = float(np.dot(final, mode)); holding_work = float(np.dot(holding, mode))
    if holding_work == 0:
        return dict(lambda_UB=None, reason='zero holding load work', physical_bound_certified=False)
    mu = (dissipation-final_work)/holding_work
    return dict(lambda_UB=1-mu if holding_work < 0 else None,
                mu_work_threshold=mu, admissible_inequality='mu >= threshold' if holding_work < 0 else 'mu <= threshold',
                internal_support_work=dissipation, final_load_work=final_work, holding_load_work=holding_work,
                physical_bound_certified=False,
                limitation='kinematically compatible Tet10 mode and GP support integral; frozen strength, quadrature model only')
