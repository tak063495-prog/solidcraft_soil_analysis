"""Independent wet-domain integration: exact quadratic roots on each ray.

Rust clips reference patches and bounds that approximation with Bernstein
coefficients. This reference instead eliminates the inner wet coordinate,
then uses adaptive 1D Gauss-Legendre between analytic topology breakpoints.
"""
from functools import lru_cache
import numpy as np


@lru_cache(None)
def rule(n):return np.polynomial.legendre.leggauss(n)


def basis(r,s,count):
    r=np.asarray(r);s=np.asarray(s);l=np.stack((1-r-s,r,s),axis=-1);dl=np.array([[-1.,-1.],[1.,0.],[0.,1.]])
    if count==3:return l,np.broadcast_to(dl,(*r.shape,3,2))
    n=np.empty((*r.shape,6));g=np.empty((*r.shape,6,2));n[...,:3]=l*(2*l-1);g[...,:3,:]=(4*l-1)[...,None]*dl
    for k,(i,j) in enumerate(((0,1),(1,2),(2,0))):
        n[...,k+3]=4*l[...,i]*l[...,j];g[...,k+3,:]=4*(l[...,i,None]*dl[j]+l[...,j,None]*dl[i])
    return n,g


def roots(a,b,c):
    scale=max(abs(a),abs(b),abs(c),np.finfo(float).tiny);a/=scale;b/=scale;c/=scale
    if abs(a)<2e-15:return [] if abs(b)<2e-15 else [-c/b]
    d=b*b-4*a*c
    if d<0:return []
    q=-.5*(b+np.copysign(np.sqrt(d),b))
    if q==0:return [-b/(2*a)]
    return [q/a,c/q]


def face_force(points,level,gamma,atol=2e-12,rtol=2e-12):
    p=np.asarray(points,dtype=float);assert p.shape in ((3,3),(6,3)) and np.isfinite(p).all() and np.isfinite(level) and gamma>0
    offset=p-p[0];z=offset[:,2]
    if len(p)==3:z=np.r_[z,(z[0]+z[1])/2,(z[1]+z[2])/2,(z[2]+z[0])/2]
    f=p[0,2]-level;a=2*(z[1]-2*z[3]);d=z[1]-a;c=2*(z[2]-2*z[5]);e=z[2]-c;b=4*z[4]-a-c-2*d-2*e
    knots=[0.,1.]
    for coefficients in ((a,d,f),(a-b+c,b-2*c+d-e,c+e+f),(b*b-4*c*a,2*b*e-4*c*d,e*e-4*c*f)):
        knots.extend(x for x in roots(*coefficients) if 0<x<1)
    knots=sorted(set(knots));evaluations=0
    inner_nodes,inner_weights=rule(4)
    def ray(r):
        nonlocal evaluations
        evaluations+=1;assert evaluations<100000,'independent wet-ray integration budget exceeded'
        length=1-r;aa=c*length*length;bb=(b*r+e)*length;cc=a*r*r+d*r+f
        stops=sorted(set([0.,1.]+[v for v in roots(aa,bb,cc) if 0<v<1]));result=np.zeros((len(p),3))
        for lo,hi in zip(stops,stops[1:]):
            mid=(lo+hi)/2
            if aa*mid*mid+bb*mid+cc>=0:continue
            t=(lo+hi)/2+(hi-lo)*inner_nodes/2;s=length*t
            n,g=basis(np.full(4,r),s,len(p));j=np.einsum('ai,qaj->qij',offset,g);normal=np.cross(j[:,:,0],j[:,:,1])
            depth=level-p[0,2]-n@offset[:,2]
            result-=np.einsum('qi,qj,q,q->ij',n,normal,gamma*depth,inner_weights*(hi-lo)*length/2)
        return result
    def integrate(lo,hi,n):
        x,w=rule(n);return sum((weight*ray((lo+hi)/2+(hi-lo)*u/2) for u,weight in zip(x,w)),np.zeros((len(p),3)))*(hi-lo)/2
    def adaptive(lo,hi,tolerance,depth):
        low=integrate(lo,hi,12);high=integrate(lo,hi,24);error=float(np.linalg.norm(high-low))
        if error<tolerance+rtol*np.linalg.norm(high):return high,error
        assert depth<24,'independent wet-ray integration did not converge'
        mid=(lo+hi)/2;left,ea=adaptive(lo,mid,tolerance/2,depth+1);right,eb=adaptive(mid,hi,tolerance/2,depth+1)
        return left+right,ea+eb
    force=np.zeros((len(p),3));error=0.
    for lo,hi in zip(knots,knots[1:]):
        if hi-lo<1e-15:continue
        value,err=adaptive(lo,hi,atol*(hi-lo),0);force+=value;error+=err
    assert np.isfinite(force).all()
    return force,{'estimated_outer_integration_error_kn':error,'ray_evaluations':evaluations,'analytic_topology_intervals':len(knots)-1}
