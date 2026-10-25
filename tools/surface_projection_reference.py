"""Independent positive-Jacobian integral: quadratic roots along rays + 1D Gauss."""
import numpy as np
from hydrostatic_reference import basis, roots, rule

def surface_area(points):
    p=np.asarray(points);previous=None
    for order in (5,10,16,24):
        x,w=rule(order);x=(x+1)/2;w=w/2;r=np.repeat(x,order);s=((1-x[:,None])*x[None,:]).ravel();weights=(w[:,None]*w[None,:]*(1-x[:,None])).ravel()
        _,g=basis(r,s,len(p));j=np.einsum('ai,qaj->qij',p-p[0],g);area=float(weights@np.linalg.norm(np.cross(j[:,:,0],j[:,:,1]),axis=1))
        if previous is not None and abs(area-previous)<1e-12*area:return area,abs(area-previous)
        previous=area
    raise AssertionError('independent surface area quadrature did not converge')

def upward_area(points):
    p=np.asarray(points);r=np.array([0.,1.,0.,.5,.5,0.]);s=np.array([0.,0.,1.,0.,.5,.5]);_,g=basis(r,s,len(p))
    j=np.einsum('ai,qaj->qij',p-p[0],g);values=np.cross(j[:,:,0],j[:,:,1])[:,2]
    bern=np.r_[values[:3],2*values[3:]-(values[[0,1,2]]+values[[1,2,0]])/2]
    if np.min(bern)>=0:return float(values[3:].sum()/6),0.
    if np.max(bern)<=0:return 0.,0.
    z=values-values[0];f=values[0];a=2*(z[1]-2*z[3]);d=z[1]-a;c=2*(z[2]-2*z[5]);e=z[2]-c;b=4*z[4]-a-c-2*d-2*e
    knots=[0.,1.]
    for coefficients in ((a,d,f),(a-b+c,b-2*c+d-e,c+e+f),(b*b-4*c*a,2*b*e-4*c*d,e*e-4*c*f)):
        knots.extend(x for x in roots(*coefficients) if 0<x<1)
    knots=sorted(set(knots));calls=0
    def ray(r):
        nonlocal calls
        calls+=1;assert calls<100000
        length=1-r;aa=c*length*length;bb=(b*r+e)*length;cc=a*r*r+d*r+f
        stops=sorted(set([0.,1.]+[v for v in roots(aa,bb,cc) if 0<v<1]));total=0.
        for lo,hi in zip(stops,stops[1:]):
            mid=(lo+hi)/2
            if aa*mid*mid+bb*mid+cc<=0:continue
            total+=length*(aa*(hi**3-lo**3)/3+bb*(hi**2-lo**2)/2+cc*(hi-lo))
        return total
    def integrate(lo,hi,n):
        x,w=rule(n);return sum(weight*ray((lo+hi)/2+(hi-lo)*v/2) for v,weight in zip(x,w))*(hi-lo)/2
    scale=max(float(np.linalg.norm(j,axis=(1,2)).max())**2,1e-300)
    def adaptive(lo,hi,depth):
        low=integrate(lo,hi,12);high=integrate(lo,hi,24);error=abs(high-low)
        if error<2e-13*scale*(hi-lo)+2e-13*abs(high):return high,error
        assert depth<26
        mid=(lo+hi)/2;l,le=adaptive(lo,mid,depth+1);r,re=adaptive(mid,hi,depth+1);return l+r,le+re
    area=error=0.
    for lo,hi in zip(knots,knots[1:]):
        if hi-lo<1e-15:continue
        value,err=adaptive(lo,hi,0);area+=value;error+=err
    return area,error
