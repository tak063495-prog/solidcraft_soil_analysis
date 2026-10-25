"""Input-only 3D cell-centred FV; numerical Newton and independent BiCGSTAB.

Unit cube, reservoir H=.5 on left z<.5, optionally .25<y<.75.
Other x faces and top are no-ponding atmosphere; base/y faces no-flow.
Homogeneous isotropic VG/Mualem; no native solutions, matrices or fluxes read.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np


def conductivity(psi):
    # Direct formula, independent of the native stable logarithmic evaluation.
    m=1-1/1.6
    saturation=(1+(1.2*np.maximum(-psi,0))**1.6)**(-m)
    return 1e-6*np.sqrt(saturation)*(1-(1-saturation**(1/m))**m)**2


def bicgstab(coefficients,columns,rhs,rtol=2e-10):
    diagonal=coefficients[:,0]
    assert np.all(diagonal>0),float(diagonal.min())
    apply=lambda x: np.sum(coefficients*x[columns],axis=1)
    x=np.zeros_like(rhs);r=rhs.copy();shadow=r.copy();p=r.copy();v=np.zeros_like(r)
    rho_old=alpha=omega=1.;target=max(1e-25,rtol*np.linalg.norm(rhs))
    for iteration in range(10000):
        rho=float(shadow@r)
        assert abs(rho)>1e-290,('linear rho breakdown',iteration)
        if iteration:p=r+(rho/rho_old)*(alpha/omega)*(p-omega*v)
        phat=p/diagonal;v=apply(phat);divisor=float(shadow@v)
        assert abs(divisor)>1e-290,('linear alpha breakdown',iteration)
        alpha=rho/divisor;s=r-alpha*v
        if np.linalg.norm(s)<=target:x+=alpha*phat;break
        shat=s/diagonal;t=apply(shat);tt=float(t@t)
        assert tt>1e-290,('linear omega breakdown',iteration)
        omega=float(t@s)/tt;x+=alpha*phat+omega*shat;r=s-omega*t
        if np.linalg.norm(r)<=target:break
        assert abs(omega)>1e-290,('linear omega zero',iteration)
        rho_old=rho
    else:raise AssertionError('linear budget exhausted')
    true=np.linalg.norm(apply(x)-rhs)
    assert true<=max(2e-24,4*target),(true,target)
    return x,iteration+1,float(true)


def solve(n,finite_width):
    assert 4<=n<=80 and n%4==0
    h=1/n;z,y,x=np.indices((n,n,n));centres=(np.arange(n)+.5)/n
    height=centres[:,None,None];head=np.broadcast_to(.5-.25*centres,(n,n,n)).copy()
    reservoir=np.broadcast_to(centres[:,None]<.5,(n,n)).copy()
    if finite_width:reservoir&=(centres[None,:]>.25)&(centres[None,:]<.75)
    def flow(H):
        k=conductivity(H-height);harm=lambda a,b:2*a*b/(a+b)
        fx=np.zeros((n,n,n+1));fy=np.zeros((n,n+1,n));fz=np.zeros((n+1,n,n))
        fx[:,:,1:-1]=-harm(k[:,:,:-1],k[:,:,1:])*np.diff(H,axis=2)/h
        fy[:,1:-1,:]=-harm(k[:,:-1,:],k[:,1:,:])*np.diff(H,axis=1)/h
        fz[1:-1,:,:]=-harm(k[:-1,:,:],k[1:,:,:])*np.diff(H,axis=0)/h
        # Outward negative flux injects reservoir water. Air never injects.
        left=2*harm(k[:,:,0],1e-6)*(H[:,:,0]-.5)/h
        left_air=np.maximum(0.,2*harm(k[:,:,0],1e-6)*(H[:,:,0]-centres[:,None])/h)
        fx[:,:,0]=-np.where(reservoir,left,left_air)
        fx[:,:,-1]=np.maximum(0.,2*harm(k[:,:,-1],1e-6)*(H[:,:,-1]-centres[:,None])/h)
        fz[-1,:,:]=np.maximum(0.,2*harm(k[-1,:,:],1e-6)*(H[-1,:,:]-1)/h)
        residual=(np.diff(fx,axis=2)+np.diff(fy,axis=1)+np.diff(fz,axis=0))*h*h
        return residual,fx,fy,fz
    indices=np.arange(n**3).reshape(n,n,n);colors=(x+2*y+3*z)%7
    columns=np.tile(indices.ravel()[:,None],(1,7));valid=np.zeros_like(columns,dtype=bool);valid[:,0]=True
    offsets=[(0,0,0),(0,0,-1),(0,0,1),(0,-1,0),(0,1,0),(-1,0,0),(1,0,0)]
    for j,(dz,dy,dx) in enumerate(offsets[1:],1):
        mask=(z+dz>=0)&(z+dz<n)&(y+dy>=0)&(y+dy<n)&(x+dx>=0)&(x+dx<n)
        valid[:,j]=mask.ravel();columns[mask.ravel(),j]=indices[(z+dz)[mask],(y+dy)[mask],(x+dx)[mask]]
    history=[];linear=[];start=time.monotonic()
    for iteration in range(50):
        value,fx,fy,fz=flow(head);norm=float(np.max(np.abs(value)));history.append(norm)
        if norm<1e-17:break
        coefficients=np.zeros_like(columns,dtype=float);delta=1e-6*(1+np.abs(head));flat_colors=colors.ravel()
        for color in range(7):
            plus=head.copy();minus=head.copy();mask=colors==color;plus[mask]+=delta[mask];minus[mask]-=delta[mask]
            difference=(flow(plus)[0]-flow(minus)[0]).ravel()
            for j in range(7):
                select=valid[:,j]&(flat_colors[columns[:,j]]==color)
                coefficients[select,j]=difference[select]/(2*delta.ravel()[columns[select,j]])
        increment,it,res=bicgstab(coefficients,columns,-value.ravel());linear.append(dict(iterations=it,true_residual_l2=res))
        for half in range(20):
            fraction=2.**-half;candidate=head+fraction*increment.reshape(head.shape)
            if np.isfinite(candidate).all() and np.max(np.abs(flow(candidate)[0]))<norm*(1-1e-4*fraction):head=candidate;break
        else:raise AssertionError(('FV line search',iteration,norm))
    else:raise AssertionError('FV Newton budget exhausted')
    # Integrate all six boundary faces, including left dry-air outflow.
    net=(fx[:,:,-1].sum()-fx[:,:,0].sum()+fz[-1,:,:].sum())*h*h
    right=fx[:,:,-1].sum()*h*h
    waterline=np.full((n,n),np.nan)
    for j in range(n):
        for i in range(n):
            psi=head[:,j,i]-centres;dry=np.flatnonzero(psi<=0)
            if not len(dry):waterline[j,i]=1.;continue
            k=int(dry[0]);waterline[j,i]=centres[0] if k==0 else centres[k-1]-psi[k-1]*h/(psi[k]-psi[k-1])
    assert abs(net)<max(1e-13,n**3*1e-17),(net,norm)
    return dict(schema='soil-fem-independent-dam-fv3d/1',cells_per_axis=n,finite_width_reservoir=finite_width,
        head_m=head.tolist(),centres_m=centres.tolist(),phreatic_height_m=waterline.tolist(),
        right_outflow_m3_s=float(right),right_seepage_area_m2=float(np.count_nonzero(fx[:,:,-1]>1e-15)*h*h),
        right_face_flux_m_s=fx[:,:,-1].tolist(),transverse_flux_rms_m_s=float(np.sqrt(np.mean(fy[:,1:-1,:]**2))),
        net_boundary_error_m3_s=float(abs(net)),maximum_cell_mass_residual_m3_s=norm,iterations=iteration,
        residual_history=history,linear_history=linear,elapsed_s=time.monotonic()-start,
        scope='Independent input-only 3D FV, harmonic face K, no-ponding atmosphere, numerical full Newton and independently coded Jacobi-preconditioned BiCGSTAB. Unit homogeneous isotropic fixture only.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--cells',type=int,required=True);p.add_argument('--finite-width',action='store_true');a=p.parse_args()
    result=solve(a.cells,a.finite_width);result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
    print({k:v for k,v in result.items() if k not in ('head_m','centres_m','phreatic_height_m','right_face_flux_m_s','linear_history')})
