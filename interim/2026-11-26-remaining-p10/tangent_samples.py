"""Sample constrained, statically condensed finite-rule tangents at the 448 anchors."""
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import cho_factor,cho_solve

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"2026-11-26-remaining-p9"))
from equilibrium import rule,P8,p6


def run(sign,out):
    ref=p6.face.ref
    g=ref.geometry(*[json.loads((ref.ROOT/"examples"/(n+".json")).read_text()) for n in
        ("joint-water-column-tet10-mesh","ground-joint-water-mc","stages-joint-water-cycle")],4,6)
    rows,w,law=rule(P8/"c02"/sign)
    nodes=np.unique(np.concatenate([ids.ravel() for ids,_,_ in rows]))
    seam=(nodes[:,None]*3+np.arange(3)).ravel()
    free=np.flatnonzero(np.isnan(g.constraint))
    inner=np.setdiff1d(free,seam)
    assert len(inner)+len(seam)==len(free)==600
    cross=g.k[np.ix_(seam,inner)]
    S=g.k[np.ix_(seam,seam)]-cross@cho_solve(cho_factor(g.k[np.ix_(inner,inner)]),cross.T)
    blocks=[]
    offset=0
    for ids,shape,_ in rows:
        # A jump shape is [-N, N]. Assemble four signed blocks of the same
        # face matrix, then sum contributions from adjacent original faces.
        pos=[(3*np.searchsorted(nodes,i)[:,None]+np.arange(3)).ravel() for i in ids]
        blocks.append((pos,shape,w[offset:offset+len(shape)]))
        offset+=len(shape)
    data=np.load(P8/"c02"/sign/"kinematics.npz")
    old=np.zeros((len(w),3))
    traction=old.copy()
    result=[]
    for (phase,t,width),u in zip(data["clock"],data["u"],strict=True):
        nod=u.reshape(-1,3)
        target=np.concatenate([N@(nod[ids[1]]-nod[ids[0]]) for ids,N,_ in rows])
        count=np.maximum(1,np.ceil(np.linalg.norm(target-old,axis=1)/width)).astype(int)
        direct,_,_=ref.project(law,old,traction,target[:,None,:],count)
        if (t*64).is_integer():
            jacobians=[]
            for h in (2e-9,1e-9):
                offsets=np.concatenate([np.eye(3),-np.eye(3)])*h
                q,_,_=ref.project(law,old,traction,target[:,None,:]+offsets,count)
                tangent=((q[:,:3]-q[:,3:])/(2*h)).transpose(0,2,1)
                at=0
                J=S.copy()
                for (a,b),N,weight in blocks:
                    T=tangent[at:at+len(N)]*weight[:,None,None]
                    B=np.stack([N.T@(T[:,i,j,None]*N) for i in range(3) for j in range(3)])
                    B=B.reshape(3,3,15,15).transpose(2,0,3,1).reshape(45,45)
                    for l,r,sgn in ((a,a,1),(a,b,-1),(b,a,-1),(b,b,1)):
                        J[np.ix_(l,r)]+=sgn*B
                    at+=len(N)
                jacobians.append(J)
            result.append(dict(phase=int(phase),t=t,smallest_singular=[float(np.linalg.svd(J,compute_uv=False)[-1]) for J in jacobians],
                perturbation_frobenius=float(np.linalg.norm(jacobians[1]-jacobians[0])),
                symmetric_part_min=float(np.linalg.eigvalsh((jacobians[-1]+jacobians[-1].T)/2)[0])))
        old,traction=target,direct[:,0]
    assert len(result)==448
    (out/f"{sign}-sampled-tangents.json").write_text(json.dumps(dict(free=len(free),interior=len(inner),seam=len(seam),
        stencils_m=[2e-9,1e-9],samples=result),indent=2)+"\n")


if __name__=="__main__":
    run(sys.argv[1],Path(sys.argv[2]))
