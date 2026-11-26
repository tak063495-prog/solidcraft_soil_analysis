"""Continuous-face coercivity matrices and a residual inclusion test on displacement balls."""
import importlib.util
import inspect
import json
import os
from math import factorial, prod
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import cho_factor, cho_solve

HERE = Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("p10_face", HERE.parent/"2026-11-26-remaining-p10/face.py")
base=importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

RADIUS = 1e-8  # Euclidean norm of all free DOFs, metres
EPS = np.finfo(float).eps

class Reduction:
    def __init__(self, refined, area):
        ref = base.prev.p6.face.ref
        inputs = [json.loads((ref.ROOT/"examples"/(s+".json")).read_text()) for s in
                  ("joint-water-column-tet10-mesh","ground-joint-water-mc","stages-joint-water-cycle")]
        self.g = ref.geometry(*inputs,4,6)
        free = np.isnan(self.g.constraint)
        self.free = free
        K = self.g.k[np.ix_(free,free)]
        pairs = np.unique(np.concatenate([v[0].T for v in refined]),axis=0)
        assert len(pairs)==25 and len(np.unique(pairs))==50
        D = np.zeros((75,len(free)))
        for i,(a,b) in enumerate(pairs):
            for j in range(3):
                D[3*i+j,3*a+j]=-1
                D[3*i+j,3*b+j]=1
        assert not D[:,~free].any()
        D=D[:,free]
        self.weight=float(os.environ.get("P11_TEST_WEIGHT","1"))
        assert self.weight in (1.,4.)
        w=np.tile([1.,1.,self.weight],len(free)//3)[free]
        K=(w[:,None]*K+K.T*w)/2
        ev,Q=np.linalg.eigh(K)
        bulk_lower=(ev[0]*(1-np.linalg.norm(np.eye(len(K))-Q@Q.T))-
                    np.linalg.norm(K-(Q*ev)@Q.T)-64*len(K)**2*EPS*np.linalg.norm(K))
        self.floors={}
        for floor in ((700.,1.) if self.weight==1 else (500.,1.)):
            # Schur complement in jump coordinates. Check inverse actions via residuals.
            B=(K+K.T)/2-floor*np.eye(len(K))
            X=cho_solve(cho_factor(B),D.T)
            H=D@X
            G=np.linalg.inv((H+H.T)/2)
            h_error=np.linalg.norm(D,2)*np.linalg.norm(B@X-D.T)/(bulk_lower-floor)
            g_error=(np.linalg.norm(G,2)*h_error+np.linalg.norm(G@H-np.eye(75)))*np.linalg.norm(G,2)
            assert bulk_lower>floor and np.linalg.norm(G,2)*h_error<.001 and g_error<.1
            self.floors[floor]=(G+G.T)/2-1.*np.eye(75)
        alpha=base.prev.ALPHA
        def multinomial(a):
            return factorial(sum(a))/prod(factorial(i) for i in a)
        # Exact simplex average of two Bernstein basis functions of degree four.
        M=np.array([[multinomial(a)*multinomial(b)/multinomial(tuple(x+y for x,y in zip(a,b)))/45
                     for a in alpha] for b in alpha])
        assert np.linalg.eigvalsh(M)[0]>0 and abs(M.sum()-1)<1e-14
        lookup={tuple(pair):i for i,pair in enumerate(pairs)}
        self.blocks=[]
        start=0
        for ids,W,_,sgn,weights,_ in refined:
            W=W.reshape(-1,15,15)
            gram=W.transpose(0,2,1)@M@W
            pos=np.array([lookup[tuple(v)] for v in ids.T])
            self.blocks.append((slice(start,start+len(W)),pos,ids,sgn,
                (gram*weights[:,None,None]).reshape(-1,225), W.mean(1)*weights[:,None],
                abs(W).max(1)*weights[:,None]))
            start+=len(W)
        assert start==len(area) and abs(area.sum()-1)<1e-12
        self.chain=True

    def evaluate(self,phase,t,u,ta,tb,gl,gh,ca,closed,all_stick,countlo,counthi,bad,crossgain,
                 local,trial_low,trial_high,kn,kt,mu,radius):
        n=len(gl)
        # Independent entries xx, xy, yy, zz. For a sliding child the disk derivative
        # is positive in the perpendicular direction. Only its radial traction couples
        # to the changing normal pressure. Enclose the uncertain direction by a cone.
        localmat=np.zeros((n,6))
        stick=all_stick & closed
        localmat[stick]=[kt,0,kt,self.weight*kn,0.,0.]
        single=(counthi==1)&(gl<=0)&~stick
        mid=(trial_high+trial_low)/2
        rad=np.linalg.norm((trial_high-trial_low)/2,axis=1)
        length=np.linalg.norm(mid,axis=1)
        e=np.divide(mid,length[:,None],out=np.tile([1.,0.],(n,1)),where=length[:,None]>1e-90)
        sin=np.minimum(1.,np.divide(rad+1e-80,np.maximum(length,1e-80))*(1+1e-12)+1e-12)
        # Lower bound on tangential block along the central perpendicular via
        # P(n) >= P(e) - sin(theta) I. Works also on the sticking branch.
        rlo=np.maximum(-gh,0)*kn*mu
        perp=kt*ca*np.minimum(1.,rlo/np.maximum(length+rad,1e-80))*(1-1e-12)
        # On an open branch the full block vanishes. If closure is uncertain,
        # rlo = 0. For closing fractions the off-diagonal derivative has norm
        # <= kt |delta_t|/|delta_g|. Young's inequality consumes the normal block.
        a=mu*kn+kt*crossgain
        b=mu*kn*sin+kt*crossgain
        radial=-perp*sin-a*(a+b)/(4*kn)
        orth=perp*(1-sin)-b*(a+b)/(4*kn)
        ex,ey=e.T
        localmat[single,0]=(radial*ex**2+orth*ey**2)[single]
        localmat[single,1]=((radial-orth)*ex*ey)[single]
        localmat[single,2]=(radial*ey**2+orth*ex**2)[single]
        multi=(counthi>1)&(gl<=0)&~stick
        if self.weight>1:
            # Normal-row test weight is applied to the complete elastic matrix above.
            # On a closed endpoint the normal diagonal is known. Retain the central
            # signed coupling and spend half the diagonal on just its uncertainty.
            sure=single & (gh<0)
            eerr=np.minimum(2.,2*sin+1e-12)
            slide=(length-rad)*(1-1e-12)>np.maximum(-gl,0)*kn*mu*(1+1e-12)
            stuck=(length+rad)*(1+1e-12)<rlo*(1-1e-12)
            alpha=np.where(slide,1.,.5)
            uncertainty=np.where(slide,eerr,.5+eerr)
            alpha=np.where(stuck|(uncertainty>=1),0.,alpha)
            uncertainty=np.where(stuck,0.,np.minimum(1.,uncertainty))
            beta=mu*kn*uncertainty+kt*crossgain
            isotropic=-perp*sin-beta**2/(2*self.weight*kn)
            localmat[sure,0]=(isotropic+perp*ey**2)[sure]
            localmat[sure,1]=(-perp*ex*ey)[sure]
            localmat[sure,2]=(isotropic+perp*ex**2)[sure]
            localmat[sure,3]=self.weight*kn/2
            localmat[sure,4]=(-mu*kn*alpha*ex/2)[sure]
            localmat[sure,5]=(-mu*kn*alpha*ey/2)[sure]
        localmat[multi]=np.column_stack((-local,0*local,-local,-local,0*local,0*local))[multi]*self.weight
        A=np.zeros((75,75))
        f=(self.g.k@u+self.g.initial-self.g.load-
           ((1-t)*(self.g.initial-self.g.load) if phase==0 else 0.)).reshape(-1,3)
        err=np.zeros_like(f)
        for sl,pos,ids,sgn,mass,mean,amax in self.blocks:
            vals=(localmat[sl].T@mass).reshape(6,15,15)
            vals[4:]*=sgn
            for val,i,j in zip(vals[[0,1,1,2,3,4,4,5,5]],(0,0,1,1,2,0,2,1,2),(0,1,0,1,2,2,0,2,1)):
                A[np.ix_(3*pos+i,3*pos+j)] += val
            zl=np.minimum(gl[sl],0)*sgn*kn
            zh=np.minimum(gh[sl],0)*sgn*kn
            lo=np.column_stack((ta[sl],np.minimum(zl,zh)))
            hi=np.column_stack((tb[sl],np.maximum(zl,zh)))
            center=mean.T@((lo+hi)/2)
            uncertainty=amax.T@((hi-lo)/2)
            f[ids[0]]-=center
            f[ids[1]]+=center
            err[ids[0]]+=uncertainty
            err[ids[1]]+=uncertainty
        residual=float(np.linalg.norm((abs(f)+err).ravel()[self.free])*(1+1e-9)+1e-8)
        lower=None
        margin=None
        if not bad.any():
            for floor,G in self.floors.items():
                # Protect assembly cancellation even when a face contributes only
                # a small net matrix. 0.1 kN/m exceeds the roundoff guard below.
                scale=max(1.,float(abs(localmat).max()))
                guard=1e-8*scale+.1
                H=G+(A+A.T)/2
                margin=float(np.linalg.eigvalsh(H)[0]-guard)
                if margin>0:
                    lower=floor/self.weight
                    break
        inclusion=bool(lower and lower*radius>residual)
        self.chain=self.chain and inclusion
        return dict(phase=int(phase),t=t,lower=lower,**({"test_weight":self.weight} if self.weight!=1. else {}),schur_margin=margin,unbounded=int(bad.sum()),
                    closed_stick=int(stick.sum()),single_switch=int(single.sum()),other_multi=int(multi.sum()),
                    residual_upper_kn=residual,displacement_radius_m=radius,inclusion=inclusion,chain=self.chain)

def run(sign,out,radius):
    # Reuse the complete P10 interval update including all substep partitions and
    # its 7-point history checks. Expand Bernstein ranges before classifications.
    src=inspect.getsource(base.run)
    def replace(a,b):
        nonlocal src
        assert src.count(a)==1
        src=src.replace(a,b)
    replace('    K=bulk()', '    reducer=Reduction(refined,area)\n    K=bulk()')
    replace('        gl,gh=lo[:,0],hi[:,0]', '''        # Each accepted endpoint is allowed its own free displacement perturbation.
        jump_radius=np.sqrt(trace2)*radius*(1+1e-12)
        change_radius=jump_radius*(1 if not history else 2)
        for col in range(4):
            delta=jump_radius if col==0 else change_radius
            lo[:,col]-=delta
            hi[:,col]+=delta
        gl,gh=lo[:,0],hi[:,0]''')
    replace('        for m in range', '        trial_low=np.zeros_like(ta)\n        trial_high=np.zeros_like(tb)\n        for m in range')
    replace('                rl=np.maximum', '''                if m==1:
                    trial_low[idx],trial_high[idx]=l,h
                rl=np.maximum''')
    replace('        history.append(row)', '''        history.append(reducer.evaluate(phase,t,u,ta,tb,gl,gh,ca,closed,all_stick,countlo,counthi,
            bad,crossgain,local,trial_low,trial_high,kn,kt,mu,radius))''')
    src=src[:src.index('    np.savez_compressed')]+'''    with (out/f"{sign}-{radius:g}.jsonl").open("x") as f:
        for row in history:
            print(json.dumps(row,allow_nan=False),file=f)
'''
    ns=dict(vars(base),Reduction=Reduction,radius=radius)
    exec(compile(src,__file__,"exec"),ns)
    ns["run"](sign,out)

if __name__=="__main__":
    run(sys.argv[1],Path(sys.argv[2]),float(sys.argv[3]))
