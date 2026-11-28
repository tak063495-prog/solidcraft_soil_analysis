"""Signed differences at identical load coordinates on each accepted half."""
import json
from pathlib import Path
import sys

import numpy as np

from alternate import LEVELS,POINTS,STAGES
from score_reload import load_rows


def score(folder):
    plan=json.loads((folder/"matched-plan.json").read_text())
    ts=sorted(set.union(*map(set,plan["grids"])))
    snapshots={}
    costs={}
    z=np.load(folder/"accepted-halves.npz")
    for n in LEVELS:
        rows=load_rows(folder/f"matched-{n}-trace.jsonl.gz")
        prior=load_rows(folder/f"{n}-trace.jsonl.gz")
        costs[n]={}
        for j,name in enumerate(STAGES):
            run=[r for r in rows if r.get("stage")==name]
            original=[r for r in prior if r.get("stage")==name]
            assert run==original[:len(run)]
            cost=len(run)+sum(r.get("stage")==j and "charge" in r for r in rows)
            assert cost<=1787
            costs[n][name]=cost
        snap={r["matched"]:r for r in rows if "matched" in r}
        assert sorted(snap)==ts
        t=z[f"recontact_{n}_progress"]
        for u,r in snap.items():
            assert r["force"]<=1e-10
            j=np.searchsorted(t,u)
            if t[j]==u:
                np.testing.assert_array_equal(z[f"recontact_{n}_slip"][j],[v[1] for v in r["points"]])
        snapshots[n]=snap
    result={}
    for i,p in enumerate(POINTS):
        tt=plan["grids"][i]
        v=np.array([[snapshots[n][u]["points"][i][1] for u in tt] for n in LEVELS])
        a,b=v[0]-v[1],v[1]-v[2]
        def at(j):
            u=tt[j]
            samples=[]
            for n in LEVELS:
                r=snapshots[n][u]
                t=z[f"recontact_{n}_progress"]
                k=int(np.searchsorted(t,u))
                left=z[f"recontact_{n}_slip"][k-1,i]
                increment=np.array(r["points"][i][1])-left
                length=float(np.linalg.norm(increment))
                samples.append(dict(n=n,increment_from_half_start=increment.tolist(),
                    direction=(increment/length).tolist() if length else None,span=t[k-1:k+1].tolist(),
                    attempt_half=z[f"recontact_{n}_attempt_half"][k].tolist(),
                    half_classes=z[f"recontact_{n}_classes"][k-1:k+1,i].tolist(),
                    class_at_probe=r["points"][i][0],slip=r["points"][i][1],force=r["force"]))
            return dict(t=u,levels=samples,coarse_medium=a[j].tolist(),medium_fine=b[j].tolist(),
                        dot=float(a[j]@b[j]))
        events={}
        for name,neg in (("dot",np.sum(a*b,axis=1)<0),("x",a[:,0]*b[:,0]<0),("y",a[:,1]*b[:,1]<0)):
            on=np.flatnonzero(neg & ~np.r_[False,neg[:-1]])
            events[name]=[dict(previous=at(k-1) if k else None,observed=at(k)) for k in on]
        result[p]=dict(rows=[at(j) for j in range(len(tt))],negative_episodes=events)
    (folder/"matched-signs.json").write_text(json.dumps(dict(costs=costs,points=result),indent=2)+"\n")


if __name__=="__main__":
    score(Path(sys.argv[1]))
