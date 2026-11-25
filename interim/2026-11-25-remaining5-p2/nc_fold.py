"""Localize the weak-mode sign change on one forward increment with fixed history."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from nc_continue import Section


def main(inputs,binary,route):
    m=Section(inputs,binary)
    report=json.loads((route/"report.json").read_bytes())
    changes=[i for i,(a,b) in enumerate(zip(report["rows"],report["rows"][1:]))
             if a["signed_soft_stiffness"]*b["signed_soft_stiffness"]<0]
    i,=changes
    # Start one full step before the lower bracket to keep a nonzero material increment.
    with gzip.open(route/"accepted.jsonl.gz","rt") as f:
        saved=[json.loads(line) for line in f]
    a,parent=saved[i-2:i]
    step=report["step_m"]
    def vector(r):
        return np.r_[np.ravel(r["displacement_m"])[m.free],r["load"]]
    t=vector(parent)-vector(a)
    t/=np.linalg.norm(t)
    samples=[]
    def solve(ds):
        z=ds*t
        counts=np.ones(len(m.initial),dtype=int)
        for it in range(m.settings["maximum_newton_iterations"]):
            ev=m.evaluate(parent["active_material_states"],z[:-1],counts)
            if ev is None:
                z=ds*t
                continue
            f,k,_,plastic=ev
            lam=parent["load"]+z[-1]
            external=m.load(lam)
            r=external[m.free]-f
            g=ds-t@z
            force=np.linalg.norm(r)/np.linalg.norm(external)
            mat=np.vstack((np.column_stack((k,m.h)),t))
            if it and force<1e-13 and abs(g)<1e-12:
                break
            rhs=np.r_[r,g]
            dx=np.linalg.solve(mat,rhs)
            fr=np.linalg.norm(k@dx[:-1]+m.h*dx[-1]-r)/np.linalg.norm(r)
            lr=np.linalg.norm(mat@dx-rhs)/np.linalg.norm(rhs)
            assert max(fr,lr)<=m.settings["relative_linear_tolerance"]
            z+=dx
        else:
            raise ArithmeticError("fold localization Newton limit")
        u,s,vh=np.linalg.svd(k)
        w=u[:,-1]*np.sign(u[:,-1]@m.h)
        v=vh[-1]*np.sign(vh[-1]@t[:-1])
        signed=w@k@v
        samples.append(dict(local_arc_m=ds,load=lam,signed_soft_stiffness=signed,
                            two_smallest_kN_m=s[-2:].tolist(),projected_h_kn=w@m.h,
                            physical_force=force,plastic_points=plastic))
        return signed,(mat,z,w,v,counts)
    lo,hi=0.5*step,2.5*step
    root=brentq(lambda ds:solve(ds)[0],lo,hi,xtol=1e-11)
    _,(mat,z,w,v,counts)=solve(root)
    certificate=samples[-1].copy()
    certificate["bordered_min_singular"]=np.linalg.svd(mat,compute_uv=False)[-1]
    coefficients=[]
    for width in (2e-7,1e-7):
        sides=[m.evaluate(parent["active_material_states"],z[:-1]+sign*width*v,counts) for sign in (-1,1)]
        assert all(side is not None for side in sides)
        coefficients.append(dict(width_m=width,plastic_counts=[side[3] for side in sides],
                                 curvature_per_m2=-(w@((sides[1][1]-sides[0][1])@v))/(2*width*(w@m.h))))
    print(json.dumps(dict(parent_step=parent["step"],parent_lambda=parent["load"],
        parent_arc_m=parent["step"]*step,samples=samples,certificate=certificate,
        curvature=coefficients),indent=2,allow_nan=False))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("inputs","binary","route"):
        p.add_argument(k,type=Path)
    args=p.parse_args()
    main(args.inputs,args.binary.resolve(),args.route)
