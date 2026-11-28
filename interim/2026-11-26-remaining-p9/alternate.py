"""Audit exact shared half endpoints and the piecewise-linear accepted-half trace."""
from bisect import bisect_left
import gzip
import json
from pathlib import Path
import sys

import numpy as np

from score_reload import load_rows

STAGES=("geostatic","shear","open","recontact")
LEVELS=(256,512,1024)
POINTS=(3,5,13)
HERE=Path(__file__).resolve().parent
P7=HERE.parent/"2026-11-26-remaining-p7/c01"


def history(folder,n):
    rows=load_rows(folder/f"{n}-trace.jsonl.gz")
    prior=load_rows(P7/f"alternate{n}/trace.jsonl.gz")
    stages={}
    for j,name in enumerate(STAGES):
        actual=[r for r in rows if r.get("stage")==name]
        before=[r for r in prior if r.get("stage")==name]
        fields=("from","to","accepted","reason")
        assert [[r[k] for k in fields] for r in before]==[[r[k] for k in fields] for r in actual]
        candidates={tuple(r["accepted_candidate"][::2]):r for r in rows if r.get("stage_index")==j}
        if not candidates: continue
        coords=[]
        slip=[]
        labels=[]
        ids=[]
        for i,r in enumerate(actual,1):
            if not r["accepted"]: continue
            row=candidates[r["from"],r["to"]]
            for side,t in enumerate(row["accepted_candidate"]):
                points=[p[side] for p in row["points"]]
                f=[p[1] for p in points]
                s=[p[0] for p in points]
                if coords and side==0:
                    assert coords[-1]==t and np.array_equal(f,slip[-1]) and labels[-1]==s
                else:
                    coords.append(t)
                    slip.append(f)
                    labels.append(s)
                    ids.append((i,side))
        assert coords[0]==0 and coords[-1]==1 and min(np.diff(coords))>0
        stages[name]=dict(t=np.array(coords),s=np.array(slip),labels=labels,ids=ids)
    return stages


def interval(h,t):
    i=max(1,bisect_left(h["t"],t))
    return i,float(h["t"][i-1]),float(h["t"][i])


def signed(v):
    a,b=v[0]-v[1],v[1]-v[2]
    return dict(coarse_medium=a.tolist(),medium_fine=b.tolist(),dot=float(a@b),
                opposite_components=((a*b)<0).tolist())


def detail(hs,t,p):
    levels=[]
    values=[]
    for h in hs:
        i,l,r=interval(h,t)
        left,right=h["s"][i-1:i+1,p]
        v=right if t==r else left if t==l else left+(t-l)/(r-l)*(right-left)
        levels.append(dict(interval=[l,r],attempt_half=h["ids"][i],classification=[h["labels"][k][p] for k in (i-1,i)],
                           slip_left=left.tolist(),slip_right=right.tolist(),
                           slip_increment=(right-left).tolist(),exact=t in (l,r)))
        values.append(v)
    return dict(t=t,levels=levels,**signed(np.array(values)))


def score(folder):
    data=[history(folder,n) for n in LEVELS]
    report={}
    raw={}
    for stage in STAGES:
        if not all(stage in h for h in data): continue
        hs=[h[stage] for h in data]
        for n,h in zip(LEVELS,hs):
            raw[f"{stage}_{n}_progress"]=h["t"]
            raw[f"{stage}_{n}_slip"]=h["s"]
            raw[f"{stage}_{n}_classes"]=np.array(h["labels"],dtype=str)
            raw[f"{stage}_{n}_attempt_half"]=h["ids"]
        all_t=np.unique(np.concatenate([h["t"] for h in hs]))
        same=sorted(set.intersection(*[set(h["t"]) for h in hs]))
        stages={}
        for p,point in enumerate(POINTS):
            def scan(grid):
                values=np.array([np.stack([np.interp(grid,h["t"],h["s"][:,p,k]) for k in range(3)],axis=-1) for h in hs])
                a,b=values[0]-values[1],values[1]-values[2]
                records={}
                for key,neg in (("dot",np.sum(a*b,axis=1)<0),("x",a[:,0]*b[:,0]<0),("y",a[:,1]*b[:,1]<0)):
                    on=np.flatnonzero(neg & ~np.r_[False,neg[:-1]])
                    def entry(j):
                        prev=detail(hs,grid[j-1],p) if j else None
                        return dict(previous=prev,observed=detail(hs,grid[j],p))
                    records[key]=dict(first=entry(on[0]) if len(on) else None,
                        endpoint_episode=entry(on[-1]) if neg[-1] else None,episodes=len(on))
                return records
            stages[point]=dict(start=detail(hs,0.,p),end=detail(hs,1.,p),
                               exact_common=scan(same),union_linear=scan(all_t))
        report[stage]=dict(shared_count=len(same),union_count=len(all_t),points=stages)
    # Identify the episode connected to the endpoint, including inherited signs.
    # Use shared actual knots for the accounting. Retain all intervening
    # half-step endpoints so a long interval cannot look like a single step.
    events={}
    for i,point in enumerate(POINTS):
        records={}
        for metric in ("dot","x","y"):
            stage="recontact"
            for name in reversed(STAGES):
                e=report[name]["points"][point]["exact_common"][metric]["endpoint_episode"]
                if e is None: break
                stage=name
                if e["previous"] is not None: break
            e=report[stage]["points"][point]["exact_common"][metric]["endpoint_episode"]
            if e is None: continue
            before=e["previous"]
            if before is None: continue
            after=e["observed"]
            levels=[]
            deltas=[]
            for n,h in zip(LEVELS,data,strict=True):
                local=h[stage]
                ids=np.flatnonzero((local["t"]>before["t"])&(local["t"]<=after["t"]))
                v=local["s"][ids,i]-local["s"][ids-1,i]
                deltas.append(local["s"][ids[-1],i]-local["s"][ids[0]-1,i])
                halves=[]
                for j,dv in zip(ids,v,strict=True):
                    halves.append(dict(attempt_half=local["ids"][j],span=local["t"][j-1:j+1].tolist(),
                        status=[local["labels"][k][i] for k in (j-1,j)],
                        delta=dv.tolist(),unit_direction=(dv/np.linalg.norm(dv)).tolist() if np.linalg.norm(dv)>0 else None))
                levels.append(dict(n=n,halves=halves,net_increment=deltas[-1].tolist()))
            records[metric]=dict(stage=stage,prior=before,first_exact=after,levels=levels,
                adjacent_increment_change=signed(np.array(deltas)))
        events[point]=records
    with gzip.open(folder/"endpoint-origins.json.gz","wt") as f:
        json.dump(events,f,allow_nan=False)
        f.write("\n")
    assert "recontact" in report
    archived=np.load(HERE.parent/"2026-11-26-remaining-p8/c01/alternate/signed.npz")["endpoint"]
    for i,h in enumerate(data):
        np.testing.assert_array_equal(h["recontact"]["s"][-1],archived[i,list(POINTS)])
    np.savez_compressed(folder/"accepted-halves.npz",**raw)
    (folder/"signs.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    score(Path(sys.argv[1]))
