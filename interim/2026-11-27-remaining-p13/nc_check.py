"""Recalculate the stopped Schur checks from the recorded native double Jacobian exactly."""
from fractions import Fraction as F
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from continue_nc import old
from plans import save


def linear(M, b):
    A=[[*row, x] for row,x in zip(M,b,strict=True)]
    n=len(A)
    for i in range(n):
        k=max(range(i,n),key=lambda j:abs(A[j][i]))
        A[i],A[k]=A[k],A[i]
        a=A[i][i]
        A[i]=[x/a for x in A[i]]
        for j in range(n):
            if j!=i:
                factor=A[j][i]
                A[j]=[x-factor*y for x,y in zip(A[j],A[i],strict=True)]
    return [row[-1] for row in A]


def inspect(work,out):
    probe=work/"probe/target/release/examples"
    model=old.nc_extended.Extended(work/"inputs",probe/"p5_homogeneous")
    result={}
    lines=[json.loads(s) for s in (out/"global.jsonl").read_text().splitlines()]
    for path in ("ordinary","predicted"):
        last=[r["frozen_prediction"] for r in lines if r["path"]==path and "frozen_prediction" in r][-1]
        end=old.read(out/(path+"-end.json.gz"))
        if last["unique_width"]:
            continue
        assert end["load"]==last["load"]
        flat=np.zeros_like(model.base)
        flat[model.free]=end["secant"][:-1]
        slope=-np.einsum("eqai,ei->eqa",model.b,flat[model.dofs]).reshape(-1,6)
        checks=[]
        for u in last["unresolved_brackets"]:
            i,w=u["point"],u["failed_midpoint"]
            strain=slope[i]*w
            needed=max(1,int(np.ceil(np.abs(strain).sum()/model.width)))
            row=dict(mode="homogeneous",law=model.laws[model.ids[i]],state=end["active_material_states"][i],
                     strain=strain.tolist(),point=i,count=1 if needed==1 else 2*needed)
            proc=subprocess.run([probe/"p5_point"],input=json.dumps(row),env=dict(os.environ,P5_POINT_TRACE="1"),
                                capture_output=True,text=True,check=True)
            for line in proc.stderr.splitlines():
                j=json.loads(line)
                if j.get("iteration")!=0:
                    continue
                J=[[F(v) for v in a] for a in j["jac_rows"]]
                rr=list(map(F,j["r"]))
                M=[a[:7] for a in J[:7]]
                v=linear(M,[-a[7] for a in J[:7]])
                y=linear(M,rr[:7])
                den=J[7][7]+sum(x*z for x,z in zip(J[7][:7],v,strict=True))
                rhs=-rr[7]+sum(x*z for x,z in zip(J[7][:7],y,strict=True))
                dx=F(j["dx"][7])
                passes=abs(rhs/den-dx)<=F(1e-14)+F(5e-9)*abs(dx)
                checks.append(dict(point=i,width=w,matrix=j["jac_rows"],rhs_vector=j["r"],native_dx=float(dx),
                    exact_schur=float(den),exact_rhs=float(rhs),exact_direction=float(rhs/den),
                    relative_difference=float(abs((rhs/den-dx)/dx)),passes=passes))
        result[path]=checks
    save(out/"direction-check.json",result)


if __name__=="__main__":
    inspect(Path(sys.argv[1]),Path(sys.argv[2]))
