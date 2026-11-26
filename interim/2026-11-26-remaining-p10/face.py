"""Directional enclosures on a twice subdivided triangle and guarded tangent tests."""
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"2026-11-26-remaining-p9"))
import face_bounds as prev
from equilibrium import rule

EPS=np.finfo(float).eps


def closed_length(a,b):
    # Monotone decreasing in both endpoint gaps.
    both=(a<=0)&(b<=0)
    cross=(a<=0)!=(b<=0)
    t=np.zeros_like(a)
    np.divide(np.maximum(-a,-b),np.abs(a-b),out=t,where=cross)
    return np.where(both,1.,t)


def interval_fraction(a,b,c,d):
    return np.maximum(0.,closed_length(b,d)-32*EPS),np.minimum(1.,closed_length(a,c)+32*EPS)


def directional(lo,hi,rl,rh):
    # For each signed component the disk projection is monotone in that
    # component. Its magnitude decreases with |perpendicular component|.
    near=np.maximum(np.maximum(lo,-hi),0.)[:,::-1]
    far=np.maximum(-lo,hi)[:,::-1]
    def edge(x,r,y):
        length=np.hypot(x,y)
        return x*np.divide(r,length,out=np.ones_like(x),where=length>r)
    low=edge(lo,np.where(lo<0,rh,rl),np.where(lo<0,near,far))
    high=edge(hi,np.where(hi>0,rh,rl),np.where(hi>0,near,far))
    pad=64*EPS*(np.maximum(-lo,hi)+rh)+1e-95
    return np.maximum(-rh,low-pad),np.minimum(rh,high+pad)


def bulk():
    ref=prev.p6.face.ref
    inputs=[json.loads((ref.ROOT/"examples"/(name+".json")).read_text()) for name in
        ("joint-water-column-tet10-mesh","ground-joint-water-mc","stages-joint-water-cycle")]
    g=ref.geometry(*inputs,4,6)
    free=np.isnan(g.constraint)
    K=g.k[np.ix_(free,free)]
    S=(K+K.T)/2
    v,Q=np.linalg.eigh(S)
    # Residual, orthogonality and product-roundoff allowance on the assembled K.
    residue=np.linalg.norm(S-(Q*v)@Q.T)
    orth=np.linalg.norm(np.eye(len(v))-Q@Q.T)
    guard=64*len(v)**2*EPS*np.linalg.norm(S)
    return dict(free=int(sum(free)),minimum=float(v[0]),lower=float(v[0]-abs(v[0])*orth-residue-guard),
                reconstruction=float(residue),roundoff_allowance=float(guard))


def run(sign,out):
    folder=prev.P8/"c02"/sign
    blocks,law=prev.geometry(folder)
    refined=[]
    for ids,_,_,areas,tiles,sgn in blocks:
        children=np.array([c for tri in tiles for c in prev.p6.split(tri)])
        W=prev.transform(children)
        # |J(x)|^2 <= 2 sum_j max_b |Bernstein(phi_j)|^2.
        trace2=2*np.sum(np.max(abs(W),axis=1)**2,axis=1)*(1+256*EPS)
        refined.append((ids,W.reshape(-1,15),np.max(abs(W).sum(axis=2)),sgn,
                        np.repeat(areas/4,4),trace2))
    area=np.concatenate([b[-2] for b in refined])
    trace2=np.concatenate([b[-1] for b in refined])
    n=len(area)
    assert n==380*16 and np.isclose(area.sum(),1.)
    data=np.load(folder/"kinematics.npz")
    clock,U=data["clock"],data["u"]
    q=np.load(folder/"slip.npz")
    qslip,weight=q["slip"],q["weight"].reshape(380,7)
    kn,kt,mu=[law[k] for k in ("normal_stiffness_kpa_m","tangential_stiffness_kpa_m","friction_coefficient")]
    active,_,law7=rule(folder)
    bary,_=prev.p6.face.triangle_rule(1)
    unit=np.array([[0.,0.],[1.,0.],[0.,1.]])
    child_index=np.full(7,-1)
    for j,tri in enumerate(c for sub in prev.p6.split(unit) for c in prev.p6.split(sub)):
        # Each of the seven symmetric points lies strictly inside a child.
        uv=np.linalg.solve((tri[1:]-tri[0]).T,(bary[:,1:]-tri[0]).T).T
        inside=(uv.min(1)>0)&(uv.sum(1)<1)
        assert np.all(child_index[inside]==-1)
        child_index[inside]=j
    assert np.all(child_index>=0)
    take=(16*np.arange(380)[:,None]+child_index).ravel()
    old7=np.zeros((2660,3))
    tau7=old7.copy()
    trace_norm=float(np.sqrt(trace2.max())*(1+32*EPS))
    trace_integral=float(area@np.sqrt(trace2)*(1+1e-10))
    K=bulk()
    ta,tb=np.zeros((n,2)),np.zeros((n,2))
    sl,su=ta.copy(),tb.copy()
    ga,gb=np.zeros(n),np.zeros(n)
    uold=np.zeros_like(U[0])
    history=[]
    anchors={k:[] for k in ("slip_error","traction_radius","tangent_radius","disk_radius","center_slip","slip_radius")}
    for (phase,t,width),u in zip(clock,U,strict=True):
        ranges=[]
        for ids,W,wn,sgn,_,_ in refined:
            a=uold.reshape(-1,3)[ids]
            b=u.reshape(-1,3)[ids]
            delta=(b[1]-a[1])-(b[0]-a[0])
            vec=np.column_stack(((b[1,:,2]-b[0,:,2])*sgn,delta[:,:2],delta[:,2]*sgn))
            B=(W@vec).reshape(-1,15,4)
            pad=160*EPS*wn*(np.max(abs(a))+np.max(abs(b)))+1e-90
            ranges.append(prev.widen(B.min(1)-pad,B.max(1)+pad))
        lo,hi=map(np.concatenate,zip(*ranges))
        gl,gh=lo[:,0],hi[:,0]
        da,db=lo[:,1:3],hi[:,1:3]
        mag=prev.norm_hi(np.maximum(-lo[:,1:],hi[:,1:]))
        small=np.linalg.norm(np.maximum(np.maximum(lo[:,1:],-hi[:,1:]),0),axis=1)*(1-32*EPS)
        countlo=np.maximum(1,np.ceil(small/width*(1-8*EPS))).astype(int)
        counthi=np.maximum(1,np.ceil(mag/width*(1+8*EPS))).astype(int)
        ca,cb=interval_fraction(ga,gb,gl,gh)
        l,h=prev.scale(da,db,ca[:,None],cb[:,None])
        sl,su=prev.widen(sl+l,su+h)
        outlo=np.full_like(ta,np.inf)
        outhi=-outlo
        all_stick=np.ones(n,dtype=bool)
        for m in range(1,int(counthi.max())+1):
            idx=np.flatnonzero((countlo<=m)&(m<=counthi))
            if not len(idx):
                continue
            x,y=ta[idx],tb[idx]
            a,b=ga[idx],gb[idx]
            oldl,oldh=a,b
            endlo,endhi=gl[idx],gh[idx]
            d,e=prev.widen(da[idx]/m,db[idx]/m)
            for j in range(1,m+1):
                if j==m:
                    c,f=endlo,endhi
                else:
                    # Convex-combination endpoints with explicit rounding margin.
                    pad=32*EPS*(np.maximum(abs(oldl),abs(oldh))+np.maximum(abs(endlo),abs(endhi)))+1e-95
                    c=(oldl*(m-j)+endlo*j)/m-pad
                    f=(oldh*(m-j)+endhi*j)/m+pad
                fl,fh=interval_fraction(a,b,c,f)
                vl,vh=prev.scale(d,e,fl[:,None],fh[:,None])
                vl,vh=prev.widen(kt*vl,kt*vh)
                l,h=prev.widen(x+vl,y+vh)
                rl=np.maximum(-f,0)*kn*mu*(1-16*EPS)
                rh=np.maximum(-c,0)*kn*mu*(1+16*EPS)+1e-95
                all_stick[idx]&=(f<0)&(b<=0)&(prev.norm_hi(np.maximum(-l,h))<rl)
                x,y=directional(l,h,rl[:,None],rh[:,None])
                x[c>0]=0
                y[c>0]=0
                a,b=c,f
            outlo[idx]=np.minimum(outlo[idx],x)
            outhi[idx]=np.maximum(outhi[idx],y)
        assert np.all(outlo<=outhi)
        ta,tb=outlo,outhi
        nod=u.reshape(-1,3)
        new7=np.concatenate([shape@(nod[ids[1]]-nod[ids[0]]) for ids,shape,_ in active])
        counts7=np.maximum(1,np.ceil(np.linalg.norm(new7-old7,axis=1)/width)).astype(int)
        direct,_,_=prev.p6.face.ref.project(law7,old7,tau7,new7[:,None,:],counts7)
        old7,tau7=new7,direct[:,0]
        assert np.all(ta[take]<=tau7[:,:2]) and np.all(tau7[:,:2]<=tb[take])
        # With fixed count, projection is 1-Lipschitz in its trial vector
        # and radius. Crossing fraction sensitivity is bounded by the inverse
        # substep gap change if both endpoint gaps can be perturbed.
        closed=(gb<0)&(gh<0)
        opened=(ga>0)&(gl>0)
        dg=np.maximum(np.maximum(lo[:,3],-hi[:,3]),0)
        crossing=~(closed|opened)
        uncert=crossing&(dg==0)
        change=(countlo!=counthi)
        dt=prev.norm_hi(np.maximum(-da,db))
        crossgain=np.divide(dt,dg,out=np.zeros_like(dg),where=crossing&(dg>0))
        local=(kn+2*kt+counthi*(mu*kn+2*kt*crossgain))*(1+128*EPS)
        # If end is strictly open, final traction is identically zero.
        # Fully sticking closed steps have positive semidefinite tangent.
        bad=(uncert|change)&(gl<=0)&~(all_stick&closed)
        penalty=np.where((gl>0)|(all_stick&closed),0,local)
        lower=K["lower"]-(area*trace2)@penalty*(1+1e-10)
        traction_gain=np.where(gl>0,0.,local)
        row=dict(phase=int(phase),t=t,counts_max=int(counthi.max()),count_cross_area=float(area@change),
            grazing_area=float(area@uncert),unbounded_tangent_area=float(area@bad),
            coercivity_lower=float(lower) if not bad.any() else None,
            inverse_candidate=float(1/lower) if not bad.any() and lower>0 else None,
            previous_traction_gain=int(np.any(gl<=0)),jump_to_traction_gain=float(traction_gain.max()) if not bad.any() else None,
            # Slip also stores the closed part of the relative translation.
            closed_translation_gain=None if uncert.any() else float(2+2*np.max(crossgain)*(1+128*EPS)),
            previous_slip_gain=1)
        # Conditional error transport for (free u, previous tau sup, slip
        # sup), in their physical units. Include the previous-u contribution
        # to assembled traction and both endpoints in the closed path length.
        inv=row["inverse_candidate"]
        gain=row["jump_to_traction_gain"]
        cg=row["closed_translation_gain"]
        if inv is not None and gain is not None and cg is not None:
            tu=row["previous_traction_gain"]
            displacement=np.array([trace_integral*gain*trace_norm,trace_integral*tu,0.])*inv
            jumps=trace_norm*(displacement+np.array([1.,0.,0.]))
            traction=gain*jumps+np.array([0.,tu,0.])
            slip=cg*jumps+(traction+np.array([0.,1.,0.]))/kt+np.array([0.,0.,1.])
            row["history_matrix"]=(np.stack((displacement,traction,slip))*(1+1e-10)).tolist()
        else:
            row["history_matrix"]=None
        history.append(row)
        if t*64==int(t*64):
            k=len(anchors["slip_error"])
            e,f=prev.widen(ta/kt,tb/kt)
            s0,s1=prev.widen(sl-f,su-e)
            original=qslip[k,:,:2]
            # Stored increments and the telescoping expression have different
            # floating summation order. Audit their scale-dependent roundoff.
            guard=512*EPS*(1+k*16)*max(1e-4,float(np.abs(original).max()))
            assert np.all(s0[take]-guard<=original) and np.all(original<=s1[take]+guard)
            center=(s0+s1)*.5
            sr=prev.norm_hi(np.maximum(s1-center,center-s0))+128*EPS*prev.norm_hi(np.maximum(abs(s1),abs(s0)))
            csum=(area[:,None]*center).reshape(380,16,2).sum(1)
            sr=(area*sr).reshape(380,16).sum(1)*(1+128*EPS)
            stored=qslip[k].reshape(380,7,3)[:,:,:2]
            ref=(weight[:,:,None]*stored).sum(1)
            err=(sr+prev.norm_hi(csum-ref)+256*EPS*(weight*np.linalg.norm(stored,axis=2)).sum(1))*(1+128*EPS)
            prlo=np.maximum(-gh,0)*kn*(1-16*EPS)
            prhi=np.maximum(-gl,0)*kn*(1+16*EPS)
            tmag=prev.norm_hi((tb-ta)/2)*(1+128*EPS)
            t3=np.hypot(tmag,(prhi-prlo)/2)*(1+128*EPS)
            for key,v in (("slip_error",err),("center_slip",csum),("slip_radius",sr),("traction_radius",(area*t3).reshape(380,16).sum(1)),
                          ("tangent_radius",(area*tmag).reshape(380,16).sum(1)),("disk_radius",(area*prhi*mu).reshape(380,16).sum(1))):
                anchors[key].append(v)
        ga,gb=gl,gh
        uold=u
    assert len(history)==4156 and len(anchors["slip_error"])==448
    np.savez_compressed(out/f"{sign}.npz",**anchors)
    with (out/f"{sign}-intervals.jsonl").open("w") as f:
        for r in history:
            print(json.dumps(r,allow_nan=False),file=f)
    (out/f"{sign}-bulk.json").write_text(json.dumps(K,indent=2)+"\n")


if __name__=="__main__":
    run(sys.argv[1],Path(sys.argv[2]))
