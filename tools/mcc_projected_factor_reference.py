"""Residual-corrected tangent-factor reuse for the independent MCC reference.

The previous LU is only a preconditioner for the CURRENT original matrix.
All force/volume/mean-pressure residual gates remain in the base reference.
Reuse requires strictly positive definite material tangents; noncoercive paths
use fresh LU. Fixed kinematic Bbar and a previously factored coercive system
then retain the same nullspace. There is no regularization or stress change.
"""
from collections import OrderedDict
from mcc_projected_sparse_reference import SparseProjectedSystem
import numpy as np
from scipy.sparse.linalg import splu


class CorrectedFactor:
    def __init__(self, owner, key, matrix, factor):
        self.owner=owner;self.key=key;self.matrix=matrix.copy();self.factor=factor

    @property
    def L(self):return self.factor.L

    @property
    def U(self):return self.factor.U

    def solve(self, rhs):
        rhs=np.asarray(rhs)
        if not np.isfinite(rhs).all():raise ValueError('non-finite reference linear RHS')
        value=self.factor.solve(rhs)
        scale=max(float(np.linalg.norm(rhs)),1e-30);previous=float('inf')
        for correction in range(7):
            residual=self.matrix@value-rhs
            error=float(np.linalg.norm(residual)/scale)
            if np.isfinite(error):
                self.owner.maximum_preconditioned_rhs_residual=max(self.owner.maximum_preconditioned_rhs_residual,error)
            else:self.owner.nonfinite_cached_residuals+=1
            if np.isfinite(error) and error<=1e-14:return value
            if not np.isfinite(error) or correction==6 or (correction>=2 and error>=previous*.95):break
            value-=self.factor.solve(residual);self.owner.cached_correction_solves+=1;previous=error
        # A failed cached iteration does not alter any equation or acceptance
        # tolerance. Refactor the current original matrix, as the base does.
        self.factor=splu(self.matrix);self.owner.actual_tangent_factorizations+=1
        self.owner.fresh_fallbacks+=1;self.owner.remember(self.key,self.factor)
        return self.factor.solve(rhs)


class ReusedFactorProjectedSystem(SparseProjectedSystem):
    def __init__(self,*args,**kwargs):
        self._factors=OrderedDict();self._rate_key=None;self._coercive=False
        self.actual_tangent_factorizations=0;self.cached_factor_hits=0
        self.cached_correction_solves=0;self.fresh_fallbacks=0;self.noncoercive_bypasses=0
        self.nonfinite_cached_residuals=0
        self.maximum_preconditioned_rhs_residual=0.
        super().__init__(*args,**kwargs)

    def remember(self,key,factor):
        self._factors[key]=factor;self._factors.move_to_end(key)
        while len(self._factors)>2:self._factors.popitem(last=False)

    def factorize(self,matrix):
        # splu canonicalizes CSC even when it is called on the fresh path.
        # Match that original-matrix representation on cached calls as well.
        matrix.sum_duplicates()
        if self._rate_key is None or not self._coercive:
            self.actual_tangent_factorizations+=1
            if not self._coercive:self.noncoercive_bypasses+=1
            return splu(matrix)
        key=(matrix.shape,self._rate_key)
        factor=self._factors.get(key)
        if factor is None:
            factor=splu(matrix);self.actual_tangent_factorizations+=1;self.remember(key,factor)
        else:
            self._factors.move_to_end(key);self.cached_factor_hits+=1
        return CorrectedFactor(self,key,matrix,factor)

    @staticmethod
    def material_is_coercive(ep):
        finite=bool(np.isfinite(ep).all())
        if finite:
            norm=np.max(np.abs(ep),axis=(1,2))
            symmetry=np.max(np.abs(ep-ep.transpose(0,2,1)),axis=(1,2))
            values=np.linalg.eigvalsh((ep+ep.transpose(0,2,1))/2)
            return bool(np.all(symmetry<=1e-13*np.maximum(norm,1.)) and
                np.all(values[:,0]>1e-10*np.maximum(values[:,-1],1e-30)))
        return False

    def linear_rate(self,ep,force_rate):
        self._coercive=self.material_is_coercive(ep)
        try:return super().linear_rate(ep,force_rate)
        finally:self._coercive=False

    def rate(self,state,force_rate,plastic):
        self._rate_key=(len(plastic),np.packbits(plastic).tobytes())
        try:return super().rate(state,force_rate,plastic)
        finally:self._rate_key=None

    def solver_diagnostics(self):
        report=super().solver_diagnostics()
        report['linear_rate_systems']=report.pop('factorizations')
        report.update(backend='original-matrix residual-corrected reused LU; fresh LU on noncoercive material tangents or failed correction',
            actual_tangent_lu_factorizations=self.actual_tangent_factorizations,
            cached_factor_hits=self.cached_factor_hits,cached_correction_solves=self.cached_correction_solves,
            fresh_fallbacks=self.fresh_fallbacks,noncoercive_bypasses=self.noncoercive_bypasses,
            nonfinite_cached_residuals=self.nonfinite_cached_residuals,
            maximum_preconditioned_rhs_residual=self.maximum_preconditioned_rhs_residual,
            cache_entries=len(self._factors),material_coercivity_relative_floor=1e-10,
            cached_rhs_relative_tolerance=1e-14,
            original_force_volume_mean_pressure_acceptance_tolerance=1e-10)
        return report
