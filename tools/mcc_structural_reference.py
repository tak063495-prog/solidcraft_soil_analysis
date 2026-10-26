"""Independent rate MCC with uncondensed frame/plate rotations and translations.

No Rust solved state or condensed member matrix generates the solution. The
projected-volume and GP-volume formulations share the analytical MCC continuum
law and event integration; their distinct volume spaces remain explicit.
"""
import mcc_projected_sparse_reference as sparse
import numpy as np
from scipy.sparse import bmat, csc_matrix, diags
from scipy.sparse.linalg import splu


class StructuralSystem(sparse.SparseProjectedSystem):
    def __init__(self, mesh, ground, active, case, *, projected=True, integration_method='midpoint'):
        if integration_method not in ('midpoint','rk4'):raise ValueError('unknown independent integration method')
        self.integration_method=integration_method
        self.projected = True
        # Parent initialization audits its own projected spatial operators.
        # GP-volume strain is an independent unknown, not a displacement map.
        super().__init__(mesh, ground, active, case)
        self.projected = projected
        self.volume_pressure_mass = self.mass
        if not projected:
            self.ti = np.arange(self.weight.size).reshape(self.weight.shape)
            ng, nf, np1 = self.weight.size, len(self.free), len(self.keys)
            self.gt_assembly = sparse.Assembly(self.local_free, self.ti, (nf, ng))
            self.tg_assembly = sparse.Assembly(self.ti, self.local_free, (ng, nf))
            self.gp_assembly = sparse.Assembly(self.ti, self.pi, (ng, np1))
            self.volume_pressure_mass = self.gp_assembly.matrix(self.n*self.weight[:, :, None])
        self.structures = None

    def advance(self,state,displacement,force_rate,count,progress=None):
        if self.integration_method=='rk4':
            from mcc_rate_rk4 import advance_rk4
            return advance_rk4(self,state,displacement,force_rate,count,progress)
        return super().advance(state,displacement,force_rate,count,progress)

    def configure_structures(self, model, history, factor_rate):
        """Fix active members/installation references for one load interval."""
        nf = len(self.free)
        transform = np.zeros((self.xyz.size, nf))
        if self.rotation is None:
            transform[self.free, np.arange(nf)] = 1.
        else:
            for column, slot in enumerate(self.free):
                node, axis = divmod(int(slot), 3)
                transform[3*node:3*node+3, column] = self.rotation[node, :, axis]
        _, binding, free, _ = model.space(history, np.zeros_like(self.xyz), transform)
        k, _, f0, *_ = model.assemble(history, 0.)
        _, _, f1, *_ = model.assemble(history, 1.)
        load_rate = (f1-f0)*factor_rate
        self.structures = (
            csc_matrix(binding.T@k@binding), csc_matrix(binding.T@k@free),
            csc_matrix(free.T@k@free), binding.T@load_rate, free.T@load_rate)
        self.retained_structural_dofs = free.shape[1]

    def tangent_blocks(self, ep):
        if self.projected:
            return super().tangent_blocks(ep)
        ep = ep.reshape(*self.weight.shape, 6, 6)
        dd = self.dd_assembly.matrix(np.einsum(
            'eqai,eqab,eqbj,eq->eij', self.rbd, ep, self.rbd, self.weight, optimize=True))
        cm = ep[:, :, :, :3].sum(axis=3)/3
        mc = ep[:, :, :3, :].sum(axis=2)/3
        dt = self.gt_assembly.matrix(np.einsum('eqai,eqa,eq->eiq', self.rbd, cm, self.weight))
        td = self.tg_assembly.matrix(np.einsum('eqa,eqai,eq->eqi', mc, self.rbd, self.weight))
        tt = diags((cm[:, :, :3].sum(axis=2)/3*self.weight).reshape(-1), format='csc')
        return dd, dt, td, tt

    def force(self, state):
        if self.projected:
            return super().force(state)
        return sparse.rate.System.force(self, state)

    def linear_rate(self, ep, force_rate):
        dd, dt, td, tt = self.tangent_blocks(ep)
        mass, div = self.volume_pressure_mass, self.reduced_divergence
        nt, nf, np1 = tt.shape[0], len(self.free), len(self.keys)
        if self.structures is None:
            raise ValueError('configure independent structures for the load interval first')
        bb, bq, qq, fb, fq = self.structures
        nr = qq.shape[0]
        matrix = bmat([[tt, -td, None, -mass],
                       [-dt, dd+bb, bq, -div.T],
                       [None, bq.T, qq, csc_matrix((nr, np1))],
                       [-mass.T, -div, csc_matrix((np1, nr)), None]], format='csc')
        rhs = np.r_[np.zeros(nt), self.reduce(force_rate)+fb, fq, np.zeros(np1)]
        factor = splu(matrix)
        solution = factor.solve(rhs)
        self.factorizations += 1
        self.maximum_matrix_nnz = max(self.maximum_matrix_nnz, matrix.nnz)
        self.maximum_factor_nnz = max(self.maximum_factor_nnz, factor.L.nnz+factor.U.nnz)
        for correction in range(4):
            theta, reduced, rotation, pressure = np.split(solution, (nt, nt+nf, nt+nf+nr))
            residual = matrix@solution-rhs
            scale = max(np.linalg.norm(rhs), 1e-30)
            error = float(np.linalg.norm(residual)/scale)
            blocks = dict(
                force=float(np.linalg.norm(residual[nt:nt+nf+nr])/scale),
                volume=float(np.linalg.norm(residual[-np1:])/max(np.linalg.norm(mass.T@theta), np.linalg.norm(div@reduced), 1e-30)),
                mean_pressure=float(np.linalg.norm(residual[:nt])/max(np.linalg.norm(tt@theta), np.linalg.norm(td@reduced), np.linalg.norm(mass@pressure), 1e-30)))
            if correction == 0:
                for key, value in blocks.items():
                    self.maximum_unrefined_block_residuals[key] = max(self.maximum_unrefined_block_residuals[key], value)
            if np.isfinite([error, *blocks.values()]).all() and max(error, *blocks.values()) <= 1e-10:
                break
            if correction == 3:
                raise ValueError(f'independent structural mixed block residual failed: {error}, {blocks}')
            solution -= factor.solve(residual)
            self.iterative_refinement_corrections += 1
        self.maximum_true_linear_residual = max(self.maximum_true_linear_residual, error)
        for key, value in blocks.items():
            self.maximum_block_residuals[key] = max(self.maximum_block_residuals[key], value)
        coordinates = np.zeros(self.xyz.size)
        coordinates[self.free] = reduced
        du = self.lift(coordinates)
        de = -np.einsum('eqai,ei->eqa', self.bd, du[self.dofs])
        volume = (np.einsum('eqa,ea->eq', self.n, theta[self.pi])
                  if self.projected else theta.reshape(self.weight.shape))
        de[:, :, :3] += volume[:, :, None]/3
        return du, de.reshape(-1, 6)

    def solver_diagnostics(self):
        return dict(super().solver_diagnostics(),
            volume_space='material-split P1' if self.projected else 'independent integration points',
            retained_structural_dofs=self.retained_structural_dofs,
            independent_structure_condensation=False)
