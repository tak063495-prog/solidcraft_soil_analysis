"""Independent checks of the structured 3D FV benchmark implementation."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from richards_dam_reference import bicgstab,conductivity,solve


class Reference(unittest.TestCase):
    def test_nonsymmetric_stencil_solver_matches_dense_original_equations(self):
        n=19;matrix=np.diag(np.linspace(3.,5.,n))+np.diag(np.full(n-1,-.7),1)+np.diag(np.full(n-1,-1.2),-1)
        columns=np.tile(np.arange(n)[:,None],(1,3));coeff=np.zeros((n,3));coeff[:,0]=matrix.diagonal()
        columns[1:,1]=np.arange(n-1);coeff[1:,1]=-1.2
        columns[:-1,2]=np.arange(1,n);coeff[:-1,2]=-.7
        rhs=np.sin(np.arange(n)+.3);value,_,_=bicgstab(coeff,columns,rhs)
        self.assertLess(np.max(np.abs(matrix@value-rhs)),1e-9)
        np.testing.assert_allclose(value,np.linalg.solve(matrix,rhs),rtol=1e-9,atol=1e-10)

    def test_full_width_3d_matches_separately_assembled_2d_dense_newton(self):
        # A separate 2D discretization uses one finite difference per DOF,
        # a dense matrix and direct LU. No 3D colouring or BiCGSTAB is reused.
        n=8;h=1/n;z=(np.arange(n)+.5)/n;H=np.broadcast_to(.5-.25*z,(n,n)).copy()
        def residual(H):
            k=conductivity(H-z[:,None]);harm=lambda a,b:2*a*b/(a+b)
            x=np.zeros((n,n+1));v=np.zeros((n+1,n))
            x[:,1:-1]=harm(k[:,:-1],k[:,1:])*(H[:,:-1]-H[:,1:])/h
            v[1:-1]=harm(k[:-1],k[1:])*(H[:-1]-H[1:])/h
            x[:,0]=np.where(z<.5,2*harm(k[:,0],1e-6)*(.5-H[:,0])/h,-np.maximum(0,2*harm(k[:,0],1e-6)*(H[:,0]-z)/h))
            x[:,-1]=np.maximum(0,2*harm(k[:,-1],1e-6)*(H[:,-1]-z)/h)
            v[-1]=np.maximum(0,2*harm(k[-1],1e-6)*(H[-1]-1)/h)
            return (np.diff(x,axis=1)+np.diff(v,axis=0))*h,x
        for _ in range(30):
            value,_=residual(H)
            if np.max(np.abs(value))<1e-17:break
            jac=np.empty((n*n,n*n))
            for i in range(n*n):
                plus=H.copy();minus=H.copy();plus.flat[i]+=1e-6;minus.flat[i]-=1e-6
                jac[:,i]=((residual(plus)[0]-residual(minus)[0])/2e-6).ravel()
            d=np.linalg.solve(jac,-value.ravel()).reshape(n,n)
            for half in range(18):
                candidate=H+d*2.**-half
                if np.max(np.abs(residual(candidate)[0]))<np.max(np.abs(value)):H=candidate;break
            else:self.fail('independent 2D line search failed')
        else:self.fail('independent 2D Newton failed')
        ref=solve(n,False);head=np.asarray(ref['head_m'])
        np.testing.assert_allclose(head,np.broadcast_to(H[:,None,:],head.shape),rtol=0,atol=2e-9)
        self.assertLess(abs(ref['right_outflow_m3_s']-residual(H)[1][:,-1].sum()*h),1e-15)
        self.assertLess(ref['transverse_flux_rms_m_s'],1e-14)


if __name__=='__main__':unittest.main()
