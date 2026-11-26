"""Compare the three actual reload histories at an identical accepted grid point."""
import gzip
import hashlib
import json
from pathlib import Path
import sys

def summarize(folder):
    paths={}
    plan=json.loads((folder/"plan.json").read_text())
    for name in ("cross","ordinary","half"):
        with gzip.open(folder/(name+".jsonl.gz"),"rt") as stream:
            rows=[json.loads(s) for s in stream if s.startswith("{")]
        path=paths[name]=[r for r in rows if r.get("stage")=="reload"]
        assert len(path)<=1787 and not any(r.get("charge") for r in rows if r.get("stage")==5)
        for r in path:
            if r["accepted"]:
                a=r["audit"]
                assert a["event"] is None and a["max_slip"]<=a["slip_budget"] and a["max_work"]<=a["work_budget"] and r["force"]<=1e-10
            else:
                assert r["parent_unchanged"]
    shared=set.intersection(*[{r["to"] for r in path if r["accepted"] and (r["to"]*1024).is_integer()} for path in paths.values()])
    def metrics(t):
        out={}
        for name,path in paths.items():
            end=len(path) if t is None else 1+next(i for i,r in enumerate(path) if r["accepted"] and r["to"]==t)
            seq=path[:end]
            a=[r for r in seq if r["accepted"]]
            out[name]=dict(accepted=len(a),rejected=end-len(a),attempts=end,progress=a[-1]["to"])
        return out
    # Reconstruct every frozen crossing action on the independently evolving path.
    import math
    increment=1/1024
    changed=0
    actions=[]
    for r in paths["cross"]:
        parent=r["from"]
        grid=(math.floor(parent*1024)+1)/1024
        if grid<=parent+16*sys.float_info.epsilon:
            grid+=1/1024
        nominal=min(1.,parent+increment,grid)
        target=nominal
        if not any(a<=parent<b for a,b in plan["bands"]):
            targets=[e["target"] for e in plan["controls"] if parent<e["span"][0] and e["span"][1]<=nominal]
            target=min([nominal,*targets])
        assert target==r["to"]
        changed+=target<nominal
        if target<nominal:
            labels=sorted({tuple(e["transition"]) for e in plan["controls"] if e["target"]==target and
                           parent<e["span"][0] and e["span"][1]<=nominal})
            actions.append(dict(parent=parent,nominal=nominal,target=target,labels=list(map(list,labels)),accepted=r["accepted"]))
        a=r["audit"]
        increment=a["proposed_increment"] if a and a["proposed_increment"] is not None else (target-parent)/2
    t=max(shared)
    # Differences are charged between common grid endpoints. Delayed path effects
    # can appear in a later cell without an action, so keep those cells as well.
    offsets={k:0 for k in paths}
    groups={}
    last=0.
    for knot in sorted(shared):
        acts=[a for a in actions if last<=a["parent"]<knot]
        kinds=sorted({"->".join(v) for a in acts for v in a["labels"]})
        group=groups.setdefault("+".join(kinds) or "no_clip",dict(cells=0,clips=0,extra_accepted=0,rejections_saved=0,extra_attempts=0))
        counts={}
        for k,seq in paths.items():
            i=offsets[k]
            while not(seq[i]["accepted"] and seq[i]["to"]==knot):
                i+=1
            segment=seq[offsets[k]:i+1]
            counts[k]=(sum(r["accepted"] for r in segment),len(segment))
            offsets[k]=i+1
        extra=counts["cross"][0]-counts["ordinary"][0]
        delta=counts["cross"][1]-counts["ordinary"][1]
        group["cells"]+=1
        group["clips"]+=len(acts)
        group["extra_accepted"]+=extra
        group["rejections_saved"]+=extra-delta
        group["extra_attempts"]+=delta
        last=knot
    return dict(actions=actions,by_controlled_cell=groups,common=metrics(t),prior_coordinate=metrics(.8466796875) if .8466796875 in shared else None,
                budget=metrics(None),shortened=changed,plan_sha256=hashlib.sha256((folder/"plan.json").read_bytes()).hexdigest())

if __name__=="__main__":
    folder=Path(sys.argv[1])
    (folder/"comparison.json").write_text(json.dumps(summarize(folder),indent=2)+"\n")
