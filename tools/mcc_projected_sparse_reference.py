"""Independent sparse rate reference for material-split P1 projected MCC.

Solves the uncondensed (compression volume, displacement, pressure) system.
No regularization, modified material law, Rust state seed or dense fallback.
The event integrator and analytic continuum MCC equations remain unchanged.
"""
import mcc_rate_fem as rate  # Set BLAS thread limits before importing NumPy.
import numpy as np
import scipy
from scipy.sparse import coo_matrix, bmat, csc_matrix
from scipy.sparse.linalg import splu
import fem_reference as ref


class Assembly:
    """Cached duplicate accumulation for element matrices, including fixed DOFs."""

    def __init__(self, rows, columns, shape):
        self.shape = shape
        local_shape = (len(rows), rows.shape[1], columns.shape[1])
        r = np.broadcast_to(rows[:, :, None], local_shape).reshape(-1)
        c = np.broadcast_to(columns[:, None, :], local_shape).reshape(-1)
        self.mask = (r >= 0) & (c >= 0)
        keys, self.inverse = np.unique(r[self.mask]*shape[1]+c[self.mask], return_inverse=True)
        self.rows, self.columns = keys//shape[1], keys % shape[1]

    def matrix(self, local):
        values = np.bincount(self.inverse, weights=local.reshape(-1)[self.mask], minlength=len(self.rows))
        return coo_matrix((values, (self.rows, self.columns)), shape=self.shape).tocsc()


class SparseProjectedSystem(rate.System):
    def __init__(self, mesh, ground, active, case=None):
        super().__init__(mesh, ground, active, case, maximum_free_dofs=100000)
        self.keys = sorted({(node, mesh['material_ids'][int(e)])
                            for e in self.element_ids for node in mesh['tetrahedra'][int(e)]})
        lookup = {key: i for i, key in enumerate(self.keys)}
        self.pi = np.asarray([[lookup[(node, mesh['material_ids'][int(e)])]
                               for node in mesh['tetrahedra'][int(e)]] for e in self.element_ids])
        q = [ref.volume(mesh, int(e)) for e in self.element_ids]
        self.n = np.stack([np.column_stack((1-v['natural'].sum(axis=1), v['natural'])) for v in q])
        reduced = np.full(self.xyz.size, -1)
        reduced[self.free] = np.arange(len(self.free))
        self.local_free = reduced[self.dofs]
        nt, nf = len(self.keys), len(self.free)
        self.tt_assembly = Assembly(self.pi, self.pi, (nt, nt))
        self.td_assembly = Assembly(self.pi, self.local_free, (nt, nf))
        self.dt_assembly = Assembly(self.local_free, self.pi, (nf, nt))
        self.dd_assembly = Assembly(self.local_free, self.local_free, (nf, nf))
        full_assembly = Assembly(self.pi, self.dofs, (nt, self.xyz.size))
        self.mass = self.tt_assembly.matrix(np.einsum('eqa,eqb,eq->eab', self.n, self.n, self.weight))
        self.mass_lu = splu(self.mass)
        self.divergence = full_assembly.matrix(np.einsum(
            'eqa,eqi,eq->eai', self.n, self.b[:, :, :3].sum(axis=2), self.weight))
        self.reduced_divergence = self.td_assembly.matrix(np.einsum(
            'eqa,eqi,eq->eai', self.n, self.rb[:, :, :3].sum(axis=2), self.weight))
        self.bd = self.b.copy()
        self.bd[:, :, :3] -= self.b[:, :, :3].sum(axis=2)[:, :, None, :]/3
        self.rbd = self.rb.copy()
        self.rbd[:, :, :3] -= self.rb[:, :, :3].sum(axis=2)[:, :, None, :]/3
        self.maximum_true_linear_residual = 0.
        self.maximum_block_residuals = dict(force=0., volume=0., mean_pressure=0.)
        self.maximum_unrefined_block_residuals = dict(force=0., volume=0., mean_pressure=0.)
        self.iterative_refinement_corrections = 0
        self.factorizations = 0
        self.maximum_matrix_nnz = 0
        self.maximum_factor_nnz = 0
        self.maximum_mass_residual = 0.
        self.checks = self.audit_operators()

    def mass_solve(self, rhs):
        value = self.mass_lu.solve(np.asarray(rhs))
        error = float(np.linalg.norm(self.mass@value-rhs)/max(np.linalg.norm(rhs), 1e-30))
        self.maximum_mass_residual = max(self.maximum_mass_residual, error)
        if error > 1e-11:
            raise ValueError('sparse projection mass residual exceeded tolerance')
        return value

    def volume(self, displacement):
        return self.mass_solve(self.divergence@displacement)

    def project_volume(self, theta):
        rhs = np.zeros(len(self.keys))
        local = np.einsum('eqa,eq,eq->ea', self.n, theta, self.weight)
        np.add.at(rhs, self.pi.reshape(-1), local.reshape(-1))
        return self.mass_solve(rhs)

    def strain(self, displacement):
        theta = self.volume(displacement)
        out = np.einsum('eqai,ei->eqa', self.bd, displacement[self.dofs])
        out[:, :, :3] += np.einsum('eqa,ea->eq', self.n, theta[self.pi])[:, :, None]/3
        return out

    def pressure(self, state):
        mean = state[:, :3].mean(axis=1).reshape(self.weight.shape)
        rhs = np.zeros(len(self.keys))
        local = np.einsum('eqa,eq,eq->ea', self.n, mean, self.weight)
        np.add.at(rhs, self.pi.reshape(-1), local.reshape(-1))
        return self.mass_solve(rhs), rhs

    def force(self, state):
        stress = state[:, :6].reshape(*self.weight.shape, 6)
        local = -np.einsum('eqai,eqa,eq->ei', self.bd, stress, self.weight, optimize=True)
        out = np.zeros(self.xyz.size)
        np.add.at(out, self.dofs.reshape(-1), local.reshape(-1))
        return out-self.divergence.T@self.pressure(state)[0]

    def tangent_blocks(self, ep):
        ep = ep.reshape(*self.weight.shape, 6, 6)
        dd = self.dd_assembly.matrix(np.einsum(
            'eqai,eqab,eqbj,eq->eij', self.rbd, ep, self.rbd, self.weight, optimize=True))
        cm = ep[:, :, :, :3].sum(axis=3)/3
        mc = ep[:, :, :3, :].sum(axis=2)/3
        mcm = cm[:, :, :3].sum(axis=2)/3
        dt = self.dt_assembly.matrix(np.einsum(
            'eqai,eqa,eqb,eq->eib', self.rbd, cm, self.n, self.weight, optimize=True))
        td = self.td_assembly.matrix(np.einsum(
            'eqa,eqb,eqbi,eq->eai', self.n, mc, self.rbd, self.weight, optimize=True))
        tt = self.tt_assembly.matrix(np.einsum(
            'eqa,eq,eqb,eq->eab', self.n, mcm, self.n, self.weight, optimize=True))
        return dd, dt, td, tt

    def factorize(self, matrix):
        """Fresh original-matrix LU; alternative references override this hook."""
        return splu(matrix)

    def linear_rate(self, ep, force_rate):
        dd, dt, td, tt = self.tangent_blocks(ep)
        m, div = self.mass, self.reduced_divergence
        nt, nf = len(self.keys), len(self.free)
        matrix = bmat([[tt, -td, -m], [-dt, dd, -div.T], [-m, -div, None]], format='csc')
        rhs = np.concatenate((np.zeros(nt), self.reduce(force_rate), np.zeros(nt)))
        factor = self.factorize(matrix)
        solution = factor.solve(rhs)
        self.factorizations += 1
        self.maximum_matrix_nnz = max(self.maximum_matrix_nnz, matrix.nnz)
        self.maximum_factor_nnz = max(self.maximum_factor_nnz, factor.L.nnz+factor.U.nnz)
        for correction in range(4):
            theta, reduced, pressure = np.split(solution, (nt, nt+nf))
            residual = matrix@solution-rhs
            error = float(np.linalg.norm(residual)/max(np.linalg.norm(rhs), 1e-30))
            blocks = dict(
                force=float(np.linalg.norm(residual[nt:nt+nf])/max(np.linalg.norm(rhs[nt:nt+nf]), 1e-30)),
                volume=float(np.linalg.norm(residual[nt+nf:])/max(np.linalg.norm(m@theta), np.linalg.norm(div@reduced), 1e-30)),
                mean_pressure=float(np.linalg.norm(residual[:nt])/max(np.linalg.norm(tt@theta), np.linalg.norm(td@reduced), np.linalg.norm(m@pressure), 1e-30)))
            if correction == 0:
                for key, value in blocks.items():
                    self.maximum_unrefined_block_residuals[key] = max(self.maximum_unrefined_block_residuals[key], value)
            finite = bool(np.isfinite([error, *blocks.values()]).all())
            if finite and error <= 1e-10 and max(blocks.values()) <= 1e-10:
                break
            if not finite or correction == 3:
                raise ValueError(f'sparse three-field true linear residual exceeded tolerance: {error}, {blocks}')
            # Reuse the same LU factors, but form each correction RHS from the
            # original unscaled equations. No altered material, load or gate.
            solution -= factor.solve(residual)
            self.iterative_refinement_corrections += 1
        self.maximum_true_linear_residual = max(self.maximum_true_linear_residual, error)
        for key, value in blocks.items():
            self.maximum_block_residuals[key] = max(self.maximum_block_residuals[key], value)
        coordinates = np.zeros(self.xyz.size)
        coordinates[self.free] = reduced
        du = self.lift(coordinates)
        de = -np.einsum('eqai,ei->eqa', self.bd, du[self.dofs])
        de[:, :, :3] += np.einsum('eqa,ea->eq', self.n, theta[self.pi])[:, :, None]/3
        return du, de.reshape(-1, 6)

    def rate(self, state, force_rate, plastic):
        self.calls += 1
        yf, d, n, dn, den = rate.terms(state, self.parameters)
        if np.any(den[plastic] <= 0):
            raise ValueError('nonpositive active MCC consistency denominator')
        ep = d.copy()
        ep[plastic] -= np.einsum('qa,qb,q->qab', dn[plastic], dn[plastic], 1/den[plastic])
        du, de = self.linear_rate(ep, force_rate)
        divisor = np.where(plastic, den, np.maximum(np.abs(den), 1.))
        potential = np.einsum('qa,qa->q', dn, de)/divisor
        flow = (potential*plastic)[:, None]*n
        dy = np.zeros_like(state)
        dy[:, :6] = np.einsum('qab,qb->qa', d, de-flow)
        dy[:, 6] = np.exp(state[:, 7])/(self.parameters[:, 1]-self.parameters[:, 2])*flow[:, :3].sum(axis=1)
        dy[:, 7] = -de[:, :3].sum(axis=1)
        dy[:, 8:] = flow
        return du, dy, potential*np.linalg.norm(n, axis=1), yf

    def audit_operators(self):
        gradient = np.array([[.1, .2, -.3], [.4, -.2, .5], [.2, -.3, .6]])
        affine = (self.xyz@gradient.T+np.array([.2, -.5, .3])).reshape(-1)
        error = float(abs(self.strain(affine)-np.einsum('eqai,ei->eqa', self.b, affine[self.dofs])).max())
        rng = np.random.default_rng(1251025)
        stress = rng.normal(size=(*self.weight.shape, 6))
        state = np.zeros((self.weight.size, 14))
        state[:, :6] = stress.reshape(-1, 6)
        virtual = rng.normal(size=self.xyz.size)
        work = float(np.einsum('eqa,eqa,eq->', self.strain(virtual), stress, self.weight))
        adjoint = abs(float(virtual@self.force(state))+work)/max(abs(work), 1.)
        if error >= 1e-11 or adjoint >= 1e-11:
            raise ValueError('sparse projected operator audit failed')
        return dict(affine_patch_maximum_strain_error=error, virtual_work_relative_error=adjoint)

    def solver_diagnostics(self):
        return dict(backend='scipy.sparse.linalg.splu', scipy=scipy.__version__, numpy=np.__version__,
                    formulation='uncondensed compression-volume/displacement/pressure',
                    free_displacement_dofs=len(self.free), primary_volume_dofs=len(self.keys),
                    factorizations=self.factorizations, maximum_matrix_nnz=self.maximum_matrix_nnz,
                    maximum_factor_nnz=self.maximum_factor_nnz,
                    maximum_true_linear_residual=self.maximum_true_linear_residual,
                    maximum_block_residuals=self.maximum_block_residuals,
                    maximum_unrefined_block_residuals=self.maximum_unrefined_block_residuals,
                    iterative_refinement_corrections=self.iterative_refinement_corrections,
                    maximum_mass_residual=self.maximum_mass_residual, dense_fallbacks=self.dense_fallbacks)
