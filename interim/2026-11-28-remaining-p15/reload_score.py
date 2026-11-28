"""Recompute costs, verify each p15 decision against its two live parents."""
import collections
import json
from pathlib import Path
import sys

P14 = Path(__file__).resolve().parent.with_name("2026-11-27-remaining-p14")
sys.path.insert(0,str(P14))
from plans import P11,P12,P13,save
from reload_score import FIXED,counts,read


def score(out):
    files = dict(ordinary=P11/"c01/reload/ordinary.jsonl.gz",p11=P11/"c01/reload/cross.jsonl.gz",
                 p12=P12/"c01/reload/entrance.jsonl.gz",p13_replay=P13/"c01/reload/trace.jsonl.gz",
                 p14=P14/"c01/reload/trace.jsonl.gz",p15=out/"trace.jsonl.gz")
    traces = {s:read(p) for s,p in files.items()}
    paths = {s:[r for r in trace if r.get("stage")=="reload" and "accepted" in r] for s,trace in traces.items()}
    assert all(len(rows)<=1787 for rows in paths.values())
    assert not any(r.get("stage")==5 and "charge" in r for r in traces["p15"])
    trials=paths["p15"]
    forecasts={r["attempt"]:r for r in traces["p15"] if "p14_parent" in r}
    parents=[]
    hits,rejected=collections.Counter(),collections.Counter()
    left,right=.9034196486885337,.9034214389463122
    for i,r in enumerate(trials,1):
        p=r["from"]
        a=r["audit"]
        if r["accepted"]:
            assert r["force"]<=1e-10 and a["event"] is None
            assert a["max_work"]<=a["work_budget"] and a["max_slip"]<=a["slip_budget"]
        else:
            assert r["parent_unchanged"]
            if p>=left and a and a["event"]:
                rejected[a["event"][1]]+=1
        if i not in forecasts:
            assert p<.8401877295048414
            continue
        n=forecasts[i]
        assert n["p14_parent"]==p
        if not parents or p!=parents[-1]:
            parents.append(p)
        prev=parents[-2] if len(parents)>1 else None
        assert n["previous"]==prev
        proposed=n["nominal"]
        if p<left<=proposed:
            proposed=left
        elif left<=p<right:
            proposed=min(proposed,p+2*(right-left))
        for pt,old,gap,zero in n["forecasts"]:
            assert old<gap<0 and (p<left or pt in (0,28,29,31))
            distance=-gap*(p-prev)/(gap-old)
            assert zero==p+distance and zero<=n["nominal"]
            proposed=min(proposed,p+max(distance,2*(right-left)))
        assert r["to"]==n["to"]==proposed
        if proposed<n["nominal"]:
            hits[n["point"]]+=1
    targets={}
    for t in FIXED:
        targets[t]={}
        for s,rows in paths.items():
            i=next(i for i,r in enumerate(rows,1) if r["to"]==t and r["accepted"])
            targets[t][s]=counts(rows[:i])
    save(out/"comparison.json",dict(common=targets,budget={s:counts(rows) for s,rows in paths.items()},
        pretrial_records=len(forecasts),controls_by_point=hits,rejected_events_after_point30=rejected,
        retrospective_comparator="p13_replay",prior_accepted_landing_table_used=False))


if __name__=="__main__":
    score(Path(sys.argv[1]))
