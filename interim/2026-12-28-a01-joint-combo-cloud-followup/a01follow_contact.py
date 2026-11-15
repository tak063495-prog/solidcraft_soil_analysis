"""Forward analytic differentiation through the existing piecewise Coulomb path.

Three target derivatives; no differencing, smoothing, stiffness floor or solver.
At exact branch boundaries this selects the closed/loading branch and reports it.
"""
import numpy as np

class D:
    __array_priority__=1000
    def __init__(self,v,d=None):
        self.v=np.asarray(v,dtype=float)
        self.d=np.zeros(self.v.shape+(3,)) if d is None else np.asarray(d,dtype=float)
    def __add__(self,b):
        b=b if isinstance(b,D) else D(b);return D(self.v+b.v,self.d+b.d)
    __radd__=__add__
    def __neg__(self):return D(-self.v,-self.d)
    def __sub__(self,b):return self+-asdual(b)
    def __rsub__(self,b):return asdual(b)+-self
    def __mul__(self,b):
        b=asdual(b);return D(self.v*b.v,self.d*b.v[...,None]+self.v[...,None]*b.d)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=asdual(b);return D(self.v/b.v,(self.d*b.v[...,None]-self.v[...,None]*b.d)/(b.v*b.v)[...,None])
    def sum(self):return D(self.v.sum(),self.d.sum(axis=0))
    def dot(self,b):return (self*b).sum()
    def norm(self):
        v=float(np.linalg.norm(self.v));return D(v,self.v@self.d/v if v else np.zeros(3))

def asdual(v):return v if isinstance(v,D) else D(v)

def point(normal,kn,kt,mu,before,qold,rest,target,count):
    n=np.asarray(normal);before=D(before);qold=D(qold);rest=D(rest)
    start=before;end=D(target,np.eye(3));path=[];boundaries=[]
    def closed(old,oldq,after):
        delta=after-old;tauold=oldq-n*oldq.dot(n);dt=delta-n*delta.dot(n)
        trial=tauold+kt*dt;gap=(after-rest).dot(n)
        normalq=kn*gap if gap.v<=0 else D(0.)
        limit=-mu*normalq;length=trial.norm()
        if float(gap.v)==0.:boundaries.append('gap=0')
        if float(length.v)==float(limit.v):boundaries.append('friction=cap')
        if length.v>limit.v and length.v>0:
            tau=trial*(limit/length);branch='slip'
        else:tau=trial;branch='stick'
        return n*normalq+tau,branch
    for k in range(1,count+1):
        after=end if k==count else start+(end-start)*(k/count)
        g0=(before-rest).dot(n);g1=(after-rest).dot(n)
        if g0.v<=0 and g1.v<=0:q,branch=closed(before,qold,after)
        elif g0.v>0 and g1.v>0:q,branch=D(np.zeros(3)),'open'
        else:
            crossing=before+(after-before)*(g0/(g0-g1));crossing=crossing-n*(crossing-rest).dot(n)
            if g0.v<=0:q,branch=D(np.zeros(3)),'open'
            else:q,branch=closed(crossing,D(np.zeros(3)),after)
        path.append(branch);before,qold=after,q
    return q.v,q.d,dict(branch=path[-1],path=path,boundaries=boundaries)

def evaluate(system,joint,target,counts):
    values=[];tangents=[];branches=[]
    for i in range(len(target)):
        q,d,branch=point(system.normal[i],system.kn[i],system.kt[i],system.mu[i],joint['jump'][i],joint['traction'][i],joint['rest'][i],target[i],int(counts[i]))
        values.append(q);tangents.append(d);branches.append(branch)
    return np.asarray(values),np.asarray(tangents),branches

def assemble(cr,d,weight,mask=None):
    if mask is not None:cr,d,weight=cr[mask],d[mask],weight[mask]
    return np.einsum('qai,qab,qbj,q->ij',cr,d,cr,weight,optimize=True)
