"""Share a corrected LU across material active sets in one fixed KKT space.

Strictly SPD material tangents keep the same nullspace for the fixed geometry,
kinematic constraints and P1 space. Plastic active-set identity is not needed
for that property. Noncoercive tangents and failed corrections still take the
original assembled/fresh path. No physical operator or gate is changed.
"""
import numpy as np
from mcc_projected_sparse_reference import SparseProjectedSystem
from mcc_projected_matrixfree_reference import MatrixFreeFactorProjectedSystem


class CoerciveSharedFactorProjectedSystem(MatrixFreeFactorProjectedSystem):
    maximum_matrix_free_corrections=20
    def __init__(self,*args,**kwargs):
        self.previous_material_active_set=None
        self.material_active_set_changes=0
        self.active_set_change_cached_hits=0
        super().__init__(*args,**kwargs)

    def rate(self,state,force_rate,plastic):
        active=np.packbits(plastic).tobytes()
        changed=self.previous_material_active_set is not None and active!=self.previous_material_active_set
        before=self.cached_factor_hits
        self._rate_key=(len(plastic),b'fixed-coercive-kinematic-space')
        try:
            result=SparseProjectedSystem.rate(self,state,force_rate,plastic)
            if changed:
                self.material_active_set_changes+=1
                self.active_set_change_cached_hits+=self.cached_factor_hits-before
            self.previous_material_active_set=active
            return result
        finally:self._rate_key=None

    def solver_diagnostics(self):
        report=super().solver_diagnostics()
        report.update(backend='fixed-space coercive material-active-set shared LU; matrix-free original action, true-residual corrections and original fresh fallback',
            material_active_set_changes=self.material_active_set_changes,
            active_set_change_cached_hits=self.active_set_change_cached_hits,
            material_active_set_is_cache_key=False,
            geometry_support_and_P1_spaces_fixed_per_system=True,
            experimental_full_nonlinear_path_acceptance=False)
        return report
