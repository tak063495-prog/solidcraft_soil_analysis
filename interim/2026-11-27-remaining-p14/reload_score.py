"""Audit recorded pretrial proposals and compare actual accepted coordinates."""
from collections import Counter
import gzip
import json
from pathlib import Path
import sys

from plans import P11, P12, P13, save

FIXED = (.900390625,.9033203125)


def read(path):
    with gzip.open(path,"rt") as f:
        return [json.loads(line) for line in f if line.startswith("{")]


def counts(rows):
    n = sum(r["accepted"] for r in rows)
    return dict(accepted=n,rejected=len(rows)-n,attempts=len(rows),
                progress=next((r["to"] for r in reversed(rows) if r["accepted"]),0.))


def run(out):
    files = dict(ordinary=P11/"c01/reload/ordinary.jsonl.gz",p11=P11/"c01/reload/cross.jsonl.gz",
                 p12=P12/"c01/reload/entrance.jsonl.gz",p13_replay=P13/"c01/reload/trace.jsonl.gz",
                 p14=out/"trace.jsonl.gz")
    paths = {}
    for name,path in files.items():
        trace = read(path)
        paths[name] = [r for r in trace if r.get("stage")=="reload" and "accepted" in r]
        assert not any(r.get("stage")==5 and "charge" in r for r in trace)
        assert len(paths[name])<=1787
    trials = paths["p14"]
    forecasts = {r["attempt"]:r for r in trace if "p14_parent" in r}
    labels = {(r["span"][0],r["span"][2]):r for r in trace if "labels" in r}
    hits, rejects = Counter(),Counter()
    switches = []
    parents = []
    for i,r in enumerate(trials,1):
        if r["from"]>=.8401877295048414:
            note = forecasts[i]
            assert note["to"]==r["to"] and note["p14_parent"]==r["from"]
            if not parents or r["from"]!=parents[-1]:
                parents.append(r["from"])
            prev = parents[-2] if len(parents)>1 else None
            assert note["previous"]==prev
            left,right=.9034196486885337,.9034214389463122
            parent=r["from"]
            proposal=note["nominal"]
            if parent<left<=proposal:
                proposal=left
            elif left<=parent<right:
                proposal=min(proposal,parent+2*(right-left))
            for point,old,gap,zero in note["forecasts"]:
                assert old<gap<0 and prev is not None
                assert zero==r["from"]+(-gap*(r["from"]-prev)/(gap-old))
                assert zero<=note["nominal"] and (parent<left or point in (29,31))
                proposal=min(proposal,parent+max(-gap*(parent-prev)/(gap-old),2*(right-left)))
            assert note["to"]==proposal
            if note["to"]<note["nominal"]:
                hits[note["point"]]+=1
        a=r["audit"]
        if r["accepted"]:
            assert a["event"] is None and a["max_slip"]<=a["slip_budget"] and a["max_work"]<=a["work_budget"]
            assert r["force"]<=1e-10
            if r["to"]>=FIXED[-1]:
                h=labels[r["from"],r["to"]]
                for j in (0,1):
                    for k in (29,30,31):
                        if h["labels"][j][k]!=h["labels"][j+1][k]:
                            switches.append(dict(point=k,attempt=i,span=h["span"][j:j+2],labels=[x[k] for x in h["labels"][j:j+2]]))
        else:
            assert r["parent_unchanged"]
            if r["from"]>=FIXED[-1] and a and a["event"]:
                rejects[a["event"][1]]+=1
    shared=set.intersection(*[{r["to"] for r in rows if r["accepted"]} for rows in paths.values()])
    common={}
    for t in sorted(set(FIXED)|{max(shared)}):
        common[t]={}
        for name,rows in paths.items():
            end=next((i for i,r in enumerate(rows,1) if r["accepted"] and r["to"]==t),None)
            common[t][name]=counts(rows[:end]) if end else None
    save(out/"comparison.json",dict(common=common,budget={k:counts(v) for k,v in paths.items()},
        after_max_p13_common={k:counts([r for r in rows if r["from"]>=FIXED[-1]]) for k,rows in paths.items()},
        fixed_first_entrance=.8401877295048414,pretrial_records=len(forecasts),controls_by_point=hits,
        rejected_event_after_common=rejects,accepted_switches=switches,
        largest_common_accepted=max(shared),prior_accepted_landing_table_used=False,
        retrospective_comparator="p13_replay"))


if __name__ == "__main__":
    run(Path(sys.argv[1]))
