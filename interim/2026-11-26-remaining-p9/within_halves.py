"""Locate sign changes of the linear accepted-half reconstruction on its own support."""
import json
from pathlib import Path
import sys

import numpy as np

from alternate import LEVELS,POINTS,STAGES,detail


def negative_ranges(a,b,c):
    # Polynomial in the dimensionless coordinate [0,1]. Scaling preserves
    # tiny physical differences; no tolerance is applied to their signs.
    scale=max(abs(a),abs(b),abs(c))
    if scale==0: return []
    a,b,c=(np.longdouble(v)/scale for v in (a,b,c))
    roots=[]
    if a==0:
        if b!=0: roots.append(-c/b)
    else:
        disc=b*b-4*a*c
        if disc>=0:
            q=-(b+np.copysign(np.sqrt(disc),b))/2
            if q!=0: roots.extend([q/a,c/q])
            else: roots.append(-b/(2*a))
    cuts=sorted({np.longdouble(0),np.longdouble(1),*(r for r in roots if 0<r<1)})
    return [(float(l),float(r)) for l,r in zip(cuts,cuts[1:])
            if (a*((l+r)/2)+b)*((l+r)/2)+c<0]


def run(folder):
    z=np.load(folder/"accepted-halves.npz")
    result={}
    for stage in STAGES:
        hs=[]
        for n in LEVELS:
            def field(name):
                return z[f"{stage}_{n}_{name}"]
            hs.append(dict(t=field("progress"),s=field("slip"),labels=field("classes").tolist(),
                           ids=field("attempt_half").tolist()))
        grid=np.unique(np.concatenate([h["t"] for h in hs]))
        result[stage]={}
        for i,p in enumerate(POINTS):
            vals=np.array([np.stack([np.interp(grid,h["t"],h["s"][:,i,k]) for k in range(2)],axis=-1) for h in hs])
            x,y=(vals[0]-vals[1]).astype(np.longdouble),(vals[1]-vals[2]).astype(np.longdouble)
            dx,dy=np.diff(x,axis=0),np.diff(y,axis=0)
            qa,qb,qc=dx*dy,dx*y[:-1]+x[:-1]*dy,x[:-1]*y[:-1]
            data={}
            for label,k in (("dot",None),("x",0),("y",1)):
                coeff=[v.sum(axis=1) if k is None else v[:,k] for v in (qa,qb,qc)]
                episodes=[]
                for j,(a,b,c) in enumerate(zip(*coeff,strict=True)):
                    for l,r in negative_ranges(a,b,c):
                        start=grid[j]+(grid[j+1]-grid[j])*l
                        end=grid[j]+(grid[j+1]-grid[j])*r
                        if episodes and start==episodes[-1][1]:
                            episodes[-1][1]=end
                        else:
                            episodes.append([start,end,j])
                def entry(e):
                    if e is None: return None
                    start,end,j=e
                    return dict(start=start,end=end,knot_pair=grid[j:j+2].tolist(),
                                halves=detail(hs,(grid[j]+grid[j+1])/2,i)["levels"])
                data[label]=dict(first=entry(episodes[0]) if episodes else None,
                    terminal=entry(episodes[-1]) if episodes and episodes[-1][1]==1 else None,
                    episodes=len(episodes))
            result[stage][p]=data
    (folder/"within-halves.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
