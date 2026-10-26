"""Matrix-free CURRENT three-field equations with a cached-LU preconditioner.

Only a strictly coercive tangent with an existing factor is attempted without
assembly. Failed correction, new active sets and noncoercive material tangents
use the original assembled equations and fresh LU. No physical term is dropped.
"""
import mcc_projected_factor_reference as reused
import numpy as np


class MatrixFreeFactorProjectedSystem(reused.ReusedFactorProjectedSystem):
    maximum_matrix_free_corrections=6
    def __init__(self,*args,**kwargs):
        self.matrix_free_operator_actions=0;self.matrix_free_assemblies_avoided=0
        self.matrix_free_corrections=0;self.matrix_free_fallbacks=0
        super().__init__(*args,**kwargs)

    def original_action(self,ep,solution):
        """Same KKT action as [[tt,-td,-M],[-dt,dd,-D.T],[-M,-D,0]]."""
        self.matrix_free_operator_actions+=1
        nt,nf=len(self.keys),len(self.free)
        theta,u,p=np.split(solution,(nt,nt+nf))
        coordinates=np.zeros(self.xyz.size);coordinates[self.free]=u
        dev=-np.einsum('eqai,ei->eqa',self.rbd,coordinates[self.dofs])
        volume=np.einsum('eqa,ea->eq',self.n,theta[self.pi])
        vol_strain=np.zeros_like(dev);vol_strain[:,:,:3]=volume[:,:,None]/3
        tensor=ep.reshape(*self.weight.shape,6,6)
        from_u=np.einsum('eqab,eqb->eqa',tensor,dev)
        from_theta=np.einsum('eqab,eqb->eqa',tensor,vol_strain)
        def mean_rhs(stress):
            local=np.einsum('eqa,eq,eq->ea',self.n,stress[:,:,:3].mean(axis=2),self.weight)
            out=np.zeros(nt);np.add.at(out,self.pi.reshape(-1),local.reshape(-1));return out
        tt_theta=mean_rhs(from_theta);minus_td_u=mean_rhs(from_u)
        mp=self.mass@p
        local=-np.einsum('eqai,eqa,eq->ei',self.rbd,from_theta+from_u,self.weight)
        force=np.zeros(self.xyz.size);np.add.at(force,self.dofs.reshape(-1),local.reshape(-1))
        force=force[self.free]-self.reduced_divergence.T@p
        mtheta=self.mass@theta;divu=self.reduced_divergence@u
        action=np.concatenate((tt_theta+minus_td_u-mp,force,-mtheta-divu))
        return action,(tt_theta,minus_td_u,mp,mtheta,divu)

    def linear_rate(self,ep,force_rate):
        nt,nf=len(self.keys),len(self.free);size=2*nt+nf
        key=((size,size),self._rate_key)
        factor=self._factors.get(key) if self._rate_key is not None else None
        if factor is None or not self.material_is_coercive(ep):return super().linear_rate(ep,force_rate)
        self._factors.move_to_end(key)
        rhs=np.concatenate((np.zeros(nt),self.reduce(force_rate),np.zeros(nt)))
        solution=factor.solve(rhs);previous=float('inf');first_blocks=None
        for correction in range(self.maximum_matrix_free_corrections+1):
            action,parts=self.original_action(ep,solution);residual=action-rhs
            tt_theta,minus_td_u,mp,mtheta,divu=parts
            error=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-30))
            blocks=dict(force=float(np.linalg.norm(residual[nt:nt+nf])/max(np.linalg.norm(rhs[nt:nt+nf]),1e-30)),
                volume=float(np.linalg.norm(residual[nt+nf:])/max(np.linalg.norm(mtheta),np.linalg.norm(divu),1e-30)),
                mean_pressure=float(np.linalg.norm(residual[:nt])/max(np.linalg.norm(tt_theta),np.linalg.norm(minus_td_u),np.linalg.norm(mp),1e-30)))
            if first_blocks is None:first_blocks=blocks.copy()
            worst=max(error,*blocks.values())
            if np.isfinite(worst) and error<=1e-14 and max(blocks.values())<=1e-13:
                self.factorizations+=1;self.cached_factor_hits+=1;self.matrix_free_assemblies_avoided+=1
                self.maximum_true_linear_residual=max(self.maximum_true_linear_residual,error)
                for name,value in blocks.items():self.maximum_block_residuals[name]=max(self.maximum_block_residuals[name],value)
                for name,value in first_blocks.items():self.maximum_unrefined_block_residuals[name]=max(self.maximum_unrefined_block_residuals[name],value)
                theta,reduced,_=np.split(solution,(nt,nt+nf))
                coordinates=np.zeros(self.xyz.size);coordinates[self.free]=reduced;du=self.lift(coordinates)
                de=-np.einsum('eqai,ei->eqa',self.bd,du[self.dofs])
                de[:,:,:3]+=np.einsum('eqa,ea->eq',self.n,theta[self.pi])[:,:,None]/3
                return du,de.reshape(-1,6)
            if not np.isfinite(worst) or correction==self.maximum_matrix_free_corrections or (correction>=2 and worst>=previous*.95):break
            solution-=factor.solve(residual);self.matrix_free_corrections+=1;previous=worst
        self.matrix_free_fallbacks+=1
        self._factors.pop(key,None) # The fallback factors the current original matrix.
        return super().linear_rate(ep,force_rate)

    def solver_diagnostics(self):
        report=super().solver_diagnostics()
        report.update(backend='matrix-free original three-field action with corrected cached LU; original assembled/fresh fallback',
            matrix_free_operator_actions=self.matrix_free_operator_actions,
            matrix_free_assemblies_avoided=self.matrix_free_assemblies_avoided,
            matrix_free_corrections=self.matrix_free_corrections,matrix_free_fallbacks=self.matrix_free_fallbacks,
            maximum_matrix_free_corrections=self.maximum_matrix_free_corrections,
            matrix_free_rhs_relative_target=1e-14,matrix_free_block_relative_target=1e-13)
        return report
