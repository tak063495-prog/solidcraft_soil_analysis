"""Range-enclose a P4 displacement trace and the history of the p8 passive law."""
from itertools import product, permutations
from math import factorial, prod
from pathlib import Path
import json
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
P8=HERE.parent/"2026-11-26-remaining-p8"
sys.path.insert(0,str(P8))
from face_pair import p6

ALPHA=[v for v in product(range(5),repeat=3) if sum(v)==4]
EPS=np.finfo(float).eps


def transform(tiles):
    # The degree-4 blossom of a product of four affine factors is their
    # symmetrised product at four vertices. Integer arithmetic until division.
    D=256
    bary=np.concatenate((1-tiles.sum(axis=2,keepdims=True),tiles),axis=2)
    v=np.rint(D*bary).astype(np.int64)
    assert np.array_equal(v/D,bary)
    w=np.empty((len(tiles),15,15))
    for col,a in enumerate(ALPHA):
        factors=[4*v[:,:,j]-k*D for j in range(3) for k in range(a[j])]
        div=D**4*prod(factorial(j) for j in a)
        for row,b in enumerate(ALPHA):
            order=set(permutations(np.repeat(range(3),b)))
            num=sum(prod(f[:,i] for f,i in zip(factors,o,strict=True)) for o in order)
            assert np.abs(num).max()<2**53 and div*len(order)<2**53
            w[:,row,col]=num/(div*len(order))
    np.testing.assert_allclose(w.sum(axis=2),1,atol=2e-15)
    return w


def geometry(folder):
    mesh,case=[json.loads((HERE.parents[1]/"examples"/(n+".json")).read_text())
               for n in ("joint-water-column-tet10-mesh","stages-joint-water-cycle")]
    keys={}
    for cell in mesh["tetrahedra"]:
        for a in product(range(5),repeat=4):
            if sum(a)==4:
                key=tuple(sorted((v,k) for v,k in zip(cell,a) if k))
                keys.setdefault(key,len(keys))
    joint,=case["joint_system"]["joints"]
    xyz=np.array(mesh["nodes"])
    sides=[mesh["face_sets"][joint[k]] for k in ("side_a","side_b")]
    blocks=[]
    faces=json.loads((folder/"tiles.json").read_text())["faces"]
    for face in faces:
        left=sides[0][face["parent"]]
        right,=[r for r in sides[1] if set(map(tuple,xyz[r]))==set(map(tuple,xyz[left]))]
        right=[next(v for v in right if np.array_equal(xyz[v],xyz[i])) for i in left]
        tri=xyz[left]
        n=np.cross(tri[1]-tri[0],tri[2]-tri[0])
        cell,=[c for c in mesh["tetrahedra"] if set(left)<=set(c)]
        fourth,=set(cell)-set(left)
        if n@(xyz[fourth]-tri[0])>0: n=-n
        assert n[0]==n[1]==0 and n[2]!=0
        ids=np.array([[keys[tuple(sorted((i,k) for i,k in zip(side,a) if k))] for a in ALPHA]
                      for side in (left,right)])
        tiles=np.array([t for leaf in face["leaves"] for t in p6.split(np.array(leaf))])
        edges=tiles[:,1:]-tiles[:,:1]
        areas=abs(edges[:,0,0]*edges[:,1,1]-edges[:,0,1]*edges[:,1,0])*abs(n[2])/2
        w=transform(tiles)
        blocks.append((ids,w.reshape(-1,15),np.max(abs(w).sum(axis=2)),areas,tiles,np.sign(n[2])))
    return blocks,joint["material"]


def widen(lo,hi):
    return np.nextafter(lo,-np.inf),np.nextafter(hi,np.inf)


def scale(lo,hi,a,b):
    return widen(np.minimum(lo*a,lo*b),np.maximum(hi*a,hi*b))


def norm_hi(x):
    return np.sqrt(np.sum(x*x,axis=-1))*(1+8*EPS)+1e-100


def run(sign,out):
    folder=P8/"c02"/sign
    blocks,law=geometry(folder)
    archive=np.load(folder/"kinematics.npz")
    clock,U=archive["clock"],archive["u"]
    areas=np.concatenate([b[3] for b in blocks])
    N=len(areas)
    assert N==380*4 and abs(areas.sum()-1)<1e-14
    n=np.repeat([b[5] for b in blocks],[len(b[3]) for b in blocks])
    kn,kt,mu=[law[k] for k in ("normal_stiffness_kpa_m","tangential_stiffness_kpa_m","friction_coefficient")]
    z=np.zeros((N,2))
    ta,tb=z.copy(),z.copy()
    sa,sb=z.copy(),z.copy()
    ga,gb=np.zeros(N),np.zeros(N)
    prev_u=np.zeros_like(U[0])
    output={k:[] for k in ("class_code","slip_center_integral","slip_radius_known","slip_radius_uncertain",
                            "traction_center_integral","traction_radius_known","traction_radius_uncertain",
                            "uncertain_area","slip_variation_uncertain","traction_variation_uncertain","slip_rule_error_upper","traction_leaf_oscillation_upper")}
    q=np.load(folder/"slip.npz")
    qs,qw=q["slip"],q["weight"].reshape(-1,7)
    stage=0
    for (phase,t,width),u in zip(clock,U,strict=True):
        ranges=[]
        for ids,w,wn,area,tiles,sgn in blocks:
            old=prev_u.reshape(-1,3)[ids]
            cur=u.reshape(-1,3)[ids]
            nod=cur[1]-cur[0]
            delta=(cur[1]-old[1])-(cur[0]-old[0])
            vectors=np.column_stack((nod[:,2]*sgn,delta))
            b=(w@vectors).reshape(-1,15,4)
            # <= 15-term dot, rational coefficient rounding and at most three
            # subtractions; absolute error budget also covers cancellation.
            pad=128*EPS*wn*(np.max(abs(cur))+np.max(abs(old)))+1e-90
            ranges.append(widen(b.min(axis=1)-pad,b.max(axis=1)+pad))
        lo,hi=[np.concatenate(v) for v in zip(*ranges)]
        gl,gh=lo[:,0],hi[:,0]
        da,db=lo[:,1:3],hi[:,1:3]
        closed=(gb<=0)&(gh<=0)
        open_=(ga>0)&(gl>0)
        ca,cb=scale(da,db,closed[:,None].astype(float),(~open_)[:,None].astype(float))
        sa,sb=widen(sa+ca,sb+cb)
        pl,ph=widen(np.maximum(-gh,0)*kn,np.maximum(-gl,0)*kn)
        rl,rh=widen(np.maximum(pl,0)*mu,ph*mu)
        rl=np.maximum(rl,0)
        e,f=widen(kt*ca,kt*cb)
        trial_lo,trial_hi=widen(ta+e,tb+f)
        minnorm=np.sqrt(np.sum(np.maximum(np.maximum(trial_lo,-trial_hi),0)**2,axis=1))*(1-8*EPS)
        maxnorm=norm_hi(np.maximum(-trial_lo,trial_hi))
        single=norm_hi(np.maximum(-lo[:,1:],hi[:,1:]))<width
        stick=single&(gh<=0)&(maxnorm<=rl)
        slide=single&(gh<=0)&(minnorm>rh)
        # 0=open, 1=stick, 2=slide, 3=unresolved. Multiple substeps
        # are enclosed by the final Coulomb disk without presuming a status.
        code=np.where(gl>0,0,np.where(stick,1,np.where(slide,2,3))).astype("u1")
        faclo=np.maximum(0,np.nextafter(np.minimum(1,rl/maxnorm),-np.inf))
        ratio=np.divide(rh,minnorm,out=np.ones_like(rh),where=minnorm>1e-200)
        fachi=np.minimum(1,np.nextafter(np.minimum(1,ratio),np.inf))
        a,b=scale(trial_lo,trial_hi,faclo[:,None],fachi[:,None])
        ta=np.where(single[:,None],np.maximum(a,-rh[:,None]),-rh[:,None])
        tb=np.where(single[:,None],np.minimum(b,rh[:,None]),rh[:,None])
        ta[code==0]=0
        tb[code==0]=0
        assert np.all(ta<=tb)
        ga,gb=gl,gh
        prev_u=u
        if t*64!=int(t*64): continue
        # Telescoping correction/kt: plastic slip = total closed tangential
        # displacement minus final tangential traction/kt, initially both zero.
        e,f=widen(ta/kt,tb/kt)
        sl,sh=widen(sa-f,sb-e)
        nz=-n
        tl=np.column_stack((ta,np.minimum(nz*pl,nz*ph)))
        th=np.column_stack((tb,np.maximum(nz*pl,nz*ph)))
        unsure=(code==3)
        output["class_code"].append(code.reshape(-1,4))
        output["uncertain_area"].append((areas*unsure).reshape(-1,4).sum(axis=1))
        for name,l,h in (("slip",sl,sh),("traction",tl,th)):
            c=(l+h)/2
            # Pad aggregate floating-point products/sums and centre rounding.
            rad=norm_hi(np.maximum(h-c,c-l))+128*EPS*norm_hi(np.maximum(abs(l),abs(h)))
            output[name+"_variation_uncertain"].append((2*rad*unsure).reshape(-1,4).max(axis=1))
            weighted=areas*rad*(1+32*EPS)
            center=(areas[:,None]*c).reshape(380,4,-1).sum(axis=1)
            known=(weighted*(~unsure)).reshape(-1,4).sum(axis=1)
            unknown=(weighted*unsure).reshape(-1,4).sum(axis=1)
            output[name+"_center_integral"].append(center)
            output[name+"_radius_known"].append(known)
            output[name+"_radius_uncertain"].append(unknown)
            if name=="slip":
                qsum=(qs[stage].reshape(380,7,3)[:,:,:2]*qw[:,:,None]).sum(axis=1)
                discrepancy=norm_hi(center-qsum)
                guard=256*EPS*np.linalg.norm((abs(qs[stage].reshape(380,7,3))*qw[:,:,None]).sum(axis=1),axis=1)
                output["slip_rule_error_upper"].append((known+unknown+discrepancy+guard)*(1+16*EPS))
            else:
                diameter=norm_hi(h.reshape(380,4,3).max(axis=1)-l.reshape(380,4,3).min(axis=1))
                output["traction_leaf_oscillation_upper"].append(areas.reshape(-1,4).sum(axis=1)*diameter*(1+32*EPS))
        stage+=1
    assert stage==448
    out.mkdir(exist_ok=True,parents=True)
    np.savez_compressed(out/(sign+".npz"),**output,area=areas.reshape(-1,4),
                        parent=np.repeat([i for i in range(len(blocks))],[len(b[3])//4 for b in blocks]),
                        subtriangles=np.concatenate([b[4] for b in blocks]).reshape(380,4,3,2))


if __name__=="__main__":
    run(sys.argv[1],Path(sys.argv[2]))
