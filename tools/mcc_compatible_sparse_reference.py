"""Independent compatible-strain MCC with audited sparse direct linear solves.

Preserves the original full B strain, continuum tangent and event integrator.
It does not use P1 pressure or volume projection and is a separate spatial
formulation for comparing projected-volume convergence against compatibility.
"""
import mcc_rate_fem as rate
import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


class SparseCompatibleSystem(rate.System):
    def __init__(self, mesh, ground, active, case=None):
        super().__init__(mesh, ground, active, case, maximum_free_dofs=100000)
        self.factorizations = 0
        self.maximum_true_linear_residual = 0.
        self.maximum_matrix_nnz = 0
        self.maximum_factor_nnz = 0

    def linear_solve(self, local_stiffness, rhs):
        n = len(self.free)
        values = np.bincount(self.inverse, weights=local_stiffness.reshape(-1)[self.mask], minlength=len(self.sparse_rows))
        matrix = coo_matrix((values, (self.sparse_rows, self.sparse_cols)), shape=(n, n)).tocsc()
        factor = splu(matrix)
        solution = factor.solve(rhs)
        error = float(np.linalg.norm(matrix@solution-rhs)/max(np.linalg.norm(rhs), 1e-30))
        self.maximum_true_linear_residual = max(self.maximum_true_linear_residual, error)
        self.factorizations += 1
        self.maximum_matrix_nnz = max(self.maximum_matrix_nnz, matrix.nnz)
        self.maximum_factor_nnz = max(self.maximum_factor_nnz, factor.L.nnz+factor.U.nnz)
        if error > 1e-10:
            raise ValueError('compatible sparse reference true linear residual exceeded tolerance')
        return solution

    def solver_diagnostics(self):
        return dict(backend='scipy.sparse.linalg.splu', scipy=scipy.__version__, numpy=np.__version__,
                    formulation='compatible full Tet10 strain; no pressure/volume projection',
                    free_displacement_dofs=len(self.free), factorizations=self.factorizations,
                    maximum_true_linear_residual=self.maximum_true_linear_residual,
                    maximum_matrix_nnz=self.maximum_matrix_nnz, maximum_factor_nnz=self.maximum_factor_nnz,
                    dense_fallbacks=self.dense_fallbacks)
