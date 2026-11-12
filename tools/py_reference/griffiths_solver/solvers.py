"""Bounded right-preconditioned GMRES; accept only original-system residuals."""
import time

import numpy as np
from scipy.sparse import csr_matrix, eye, tril, triu
from scipy.sparse.linalg import (LinearOperator, gmres, spilu, splu,
                                 spsolve_triangular, onenormest, eigsh, eigs)

LINEAR_LIMIT = 1e-10  # Existing staged solver linear gate, never the force gate.
FORCE_LIMIT = 1e-8
METHODS = ('none_raw', 'none', 'ilu0', 'ilut', 'node_block', 'component_ilut', 'direct')


def original_residuals(a, b, x, external_force_norm):
    residual = float(np.linalg.norm(a @ x-b))
    norm_b = float(np.linalg.norm(b))
    return dict(absolute_residual=residual,
                true_linear_relative_residual=residual/norm_b if norm_b else residual,
                force_normalized_linear_residual=residual/max(external_force_norm, 1.))


def ilu0(a):
    """Natural-order no-fill/no-pivot ILU(0), including structural zeros.

    Mirrors src/linalg.rs pattern semantics. SuperLU spilu(fill_factor=1) is
    deliberately not called ILU(0): pivoting, ordering and dropping differ.
    """
    a = csr_matrix(a, copy=True)
    a.sort_indices()
    rows = [dict(zip(a.indices[a.indptr[i]:a.indptr[i+1]],
                     a.data[a.indptr[i]:a.indptr[i+1]], strict=True)) for i in range(a.shape[0])]
    for i, row in enumerate(rows):
        for j in sorted(k for k in row if k < i):
            pivot = rows[j].get(j, 0.)
            if not np.isfinite(pivot) or abs(pivot) < 1e-14:
                raise ArithmeticError(f'ILU(0) pivot failure at {j}')
            multiplier = row[j]/pivot
            row[j] = multiplier
            for k, value in rows[j].items():
                if k > j and k in row:
                    row[k] -= multiplier*value
        if not np.isfinite(list(row.values())).all() or abs(row.get(i, 0.)) < 1e-14:
            raise ArithmeticError(f'ILU(0) pivot failure at {i}')
    values = np.concatenate([[row[k] for k in a.indices[a.indptr[i]:a.indptr[i+1]]]
                             for i, row in enumerate(rows)])
    factor = csr_matrix((values, a.indices, a.indptr), shape=a.shape)
    lower, upper = tril(factor, -1).tocsr()+eye(a.shape[0], format='csr'), triu(factor).tocsr()
    def solve(v):
        return spsolve_triangular(upper, spsolve_triangular(lower, v, lower=True, unit_diagonal=True), lower=False)
    return solve, int(lower.nnz+upper.nnz)


def preconditioner(a, free, method):
    if method in ('none_raw', 'none'):
        return lambda x: x, 0
    if method == 'ilu0':
        return ilu0(a)
    if method == 'ilut':
        f = spilu(a.tocsc(), drop_tol=1e-4, fill_factor=10., permc_spec='COLAMD')
        return f.solve, int(f.L.nnz+f.U.nnz)
    if method == 'node_block':
        # Disjoint physical-node blocks retain ux/uz (or ux/uy/uz) coupling.
        _, groups = np.unique(np.asarray(free)//3, return_inverse=True)
        blocks = [np.flatnonzero(groups == i) for i in range(groups.max()+1)]
        inverses = [np.linalg.inv(a[idx][:, idx].toarray()) for idx in blocks]
        def solve(v):
            out = np.zeros_like(v)
            for idx, inv in zip(blocks, inverses, strict=True):
                out[idx] = inv @ v[idx]
            return out
        return solve, int(sum(len(idx)**2 for idx in blocks))
    if method == 'component_ilut':
        blocks = [np.flatnonzero(np.asarray(free) % 3 == i) for i in range(3)]
        blocks = [idx for idx in blocks if len(idx)]
        factors = [spilu(a[idx][:, idx].tocsc(), drop_tol=1e-4, fill_factor=10., permc_spec='COLAMD') for idx in blocks]
        def solve(v):
            out = np.zeros_like(v)
            for idx, f in zip(blocks, factors, strict=True):
                out[idx] = f.solve(v[idx])
            return out
        return solve, int(sum(f.L.nnz+f.U.nnz for f in factors))
    raise ValueError('unknown preconditioner')


def compare(a, b, free, methods=METHODS, maximum_iterations=3000, restart=60, external_force_norm=None):
    if maximum_iterations < 1 or restart < 1 or maximum_iterations > 3000 or restart > 500:
        raise ValueError('iteration budget must be 1..3000, restart 1..500')
    a, b = csr_matrix(a), np.asarray(b, float)
    free = np.asarray(free)
    if (a.shape != (len(b), len(b)) or len(free) != len(b) or not len(b)
            or not np.isfinite(a.data).all() or not np.isfinite(b).all()):
        raise ValueError('invalid system')
    if set(methods)-set(METHODS):
        raise ValueError('unknown solver method')
    force_norm = float(np.linalg.norm(b)) if external_force_norm is None else float(external_force_norm)
    if not np.isfinite(force_norm) or force_norm < 0:
        raise ValueError('invalid external force norm')
    maxima = np.asarray(abs(a).max(axis=1).toarray()).ravel()
    if np.any(maxima == 0):
        raise ValueError('zero matrix row')
    scaled = a.copy()
    scaled.data *= np.repeat(1/maxima, np.diff(scaled.indptr))
    scaled_b = b/maxima
    output = []
    for method in methods:
        row = dict(method=method, iterations=0, status='error', accepted=False,
                   linear_limit=LINEAR_LIMIT, force_limit=FORCE_LIMIT,
                   external_force_norm=force_norm, restart=restart, maximum_iterations=maximum_iterations,
                   row_scaling='none' if method in ('none_raw', 'direct') else 'row_max',
                   krylov_relative_residual_history=[])
        started = time.perf_counter()
        try:
            if not np.any(b):
                x = np.zeros_like(b)
                row.update(factor_nonzeros=0, solver_info=0, setup_seconds=0.)
            elif method == 'direct':
                factor = splu(a.tocsc())
                row.update(factor_nonzeros=int(factor.L.nnz+factor.U.nnz), setup_seconds=time.perf_counter()-started)
                x = factor.solve(b)
                row['solver_info'] = 0
            else:
                work_a, work_b = (a, b) if method == 'none_raw' else (scaled, scaled_b)
                apply, nnz = preconditioner(work_a, free, method)
                row.update(factor_nonzeros=nnz, setup_seconds=time.perf_counter()-started)
                operator = LinearOperator(a.shape, matvec=lambda v: work_a @ apply(v), dtype=float)
                history = row['krylov_relative_residual_history']
                # ||r|| <= max(row_max) ||r_scaled||. This does not relax the
                # original linear tolerance when the row norms differ.
                scale_bound = 1. if method == 'none_raw' else maxima.max()
                scaled_rtol = LINEAR_LIMIT*np.linalg.norm(b)/(scale_bound*np.linalg.norm(work_b))
                y, info = gmres(operator, work_b, rtol=scaled_rtol, atol=0., restart=restart,
                                maxiter=maximum_iterations, callback=lambda v: history.append(float(v)),
                                callback_type='legacy')  # legacy maxiter counts INNER iterations
                x = apply(y)
                row.update(iterations=len(history), solver_info=int(info))
            row.update(original_residuals(a, b, x, force_norm))
            row['accepted'] = bool(np.isfinite(x).all() and row['solver_info'] == 0
                                   and row['true_linear_relative_residual'] <= LINEAR_LIMIT
                                   and row['force_normalized_linear_residual'] <= FORCE_LIMIT)
            row['status'] = 'converged' if row['accepted'] else ('iteration_limit' if row['solver_info'] > 0 else 'residual_failed')
        except (RuntimeError, ArithmeticError, ValueError, np.linalg.LinAlgError) as error:
            row['error'] = str(error)
        row['elapsed_seconds'] = time.perf_counter()-started
        output.append(row)
    return output


def estimates(a):
    """1-norm estimate and sparse extreme eigenvalue samples, never exact kappa2.

    Non-normal matrices cannot be diagnosed from eigenvalue ratios alone. LU
    failure/ARPACK nonconvergence remain explicit diagnostic failures.
    """
    a = csr_matrix(a)
    result = dict(dimension=a.shape[0], stored_nonzeros=a.nnz,
                  relative_asymmetry=float(np.linalg.norm((a-a.T).data)/max(np.linalg.norm(a.data), 1e-300)))
    try:
        lu = splu(a.tocsc())
        inverse = LinearOperator(a.shape, matvec=lu.solve, rmatvec=lambda x: lu.solve(x, trans='T'), dtype=float)
        result['condition_1_estimate'] = float(onenormest(a)*onenormest(inverse))
        start = np.linspace(1., 2., a.shape[0])
        if a.shape[0] <= 3:
            values = np.linalg.eigvals(a.toarray())
        elif result['relative_asymmetry'] < 1e-12:
            values = np.r_[eigsh(a, k=1, which='LM', v0=start, return_eigenvectors=False, maxiter=300),
                           eigsh(a, k=1, sigma=0, OPinv=inverse, v0=start, return_eigenvectors=False, maxiter=300)]
        else:
            values = np.r_[eigs(a, k=1, which='LM', v0=start, return_eigenvectors=False, maxiter=300),
                           eigs(a, k=1, sigma=0, OPinv=inverse, v0=start, return_eigenvectors=False, maxiter=300)]
        result['eigenvalue_samples_real_imag'] = [[float(v.real), float(v.imag)] for v in values]
        result['sampled_eigenvalue_magnitude_ratio'] = float(max(abs(values))/min(abs(values)))
    except (RuntimeError, ValueError, ArithmeticError) as error:
        result['estimate_error'] = str(error)
    return result
