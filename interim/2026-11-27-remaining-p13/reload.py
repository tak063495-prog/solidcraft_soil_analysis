"""Account for post-entrance regrowth and localize point 30 on accepted half steps."""
import gzip
import json
from pathlib import Path
import sys
from plans import P11, P12, save

FIXED=(.900390625,.9033203125)
BASE=P12/"c01/reload"


def read(path):
    with gzip.open(path,"rt") as f:
        return [json.loads(s) for s in f if s.startswith("{")]


def attempts(rows):
    return [r for r in rows if r.get("stage")=="reload" and "accepted" in r]


def freeze(out):
    rows=attempts(read(BASE/"entrance.jsonl.gz"))
    actions=json.loads((BASE/"comparison.json").read_text())["actions"]
    first=actions[0][2]
    pending=[]
    groups=[]
    for r in rows:
        if r["from"]<first:
            continue
        if not r["accepted"]:
            pending.append(r)
        elif pending:
            assert all(x["from"]==r["from"] and x["to"]>r["to"] for x in pending)
            groups.append(dict(parent=r["from"],accepted=r["to"],rejected=[x["to"] for x in pending]))
            pending=[]
    targets=sorted({a[2] for a in actions}|{g["accepted"] for g in groups})
    save(out/"plan.json",dict(bands=[],controls=[dict(span=[t,t],target=t) for t in targets],
        source="p12 accepted landing after one or more post-entrance rejections, plus executed p12 entrances",
        entrance=first,groups=groups))


def counts(rows):
    return dict(accepted=sum(r["accepted"] for r in rows),rejected=sum(not r["accepted"] for r in rows),
                attempts=len(rows))


def score(out):
    identity=out/"prefix-identity.json"
    audit=json.loads(identity.read_text())
    old=json.loads((BASE/"prefix-identity.json").read_text())
    assert audit["common_canonical_sha256"]==old["common_canonical_sha256"]
    if "p13_checkpoint_sha256" not in audit:
        audit["p13_checkpoint_sha256"]=audit.pop("p12_checkpoint_sha256")
        audit["p12_checkpoint_sha256"]=old["p12_checkpoint_sha256"]
        save(identity,audit)
    opt=read(out/"trace.jsonl.gz")
    paths={k:attempts(read(p)) for k,p in (("ordinary",P11/"c01/reload/ordinary.jsonl.gz"),
        ("p11",P11/"c01/reload/cross.jsonl.gz"),("p12",BASE/"entrance.jsonl.gz"))}
    paths["p13"]=attempts(opt)
    plan=json.loads((out/"plan.json").read_text())
    prior=[x for x in paths["p12"] if x["accepted"]]
    current=[x for x in paths["p13"] if x["accepted"]]
    assert len(current)>=len(prior) and prior==current[:len(prior)]
    assert all(len(rows)<=1787 for rows in paths.values())
    assert not any(r.get("stage")==5 and "charge" in r for r in opt)
    for r in paths["p13"]:
        if r["accepted"]:
            a=r["audit"]
            assert a["event"] is None and a["max_slip"]<=a["slip_budget"] and a["max_work"]<=a["work_budget"] and r["force"]<=1e-10
        else:
            assert r["parent_unchanged"]
    half={(r["span"][0],r["span"][2]):r for r in opt if "span" in r and "labels" in r}
    transitions=[]
    tail=[]
    for idx,r in enumerate(paths["p13"],1):
        if r["from"]<FIXED[-1]:
            continue
        item=dict(attempt=idx,**r)
        if r["accepted"]:
            h=half[r["from"],r["to"]]
            labels=[v[30] for v in h["labels"]]
            item["labels_30"]=labels
            for j in (0,1):
                if labels[j]!=labels[j+1]:
                    transitions.append(dict(attempt=idx,half=j+1,span=h["span"][j:j+2],labels=labels[j:j+2]))
        tail.append(item)
    save(out/"tail.json",tail)
    common={}
    for t in FIXED:
        common[t]={}
        for label,rows in paths.items():
            end=next(i+1 for i,r in enumerate(rows) if r["accepted"] and r["to"]==t)
            common[t][label]=counts(rows[:end])
    groups=plan["groups"]
    budgets={}
    for label,rows in paths.items():
        last=next(r["to"] for r in reversed(rows) if r["accepted"])
        budgets[label]=dict(**counts(rows),progress=last,
            after_common=counts([r for r in rows if r["from"]>=FIXED[-1]]))
    # Group by the accepted parent's load. Each rejection has exactly one location;
    # flags describing its failures (slip, work, point switch) may overlap.
    ends=[plan["entrance"],max(a[2] for a in json.loads((BASE/"comparison.json").read_text())["actions"]),*FIXED]
    locations={}
    for label,rows in paths.items():
        prev_width=None
        bins=[dict(accepted=0,rejected=0,grew_and_rejected=0,switch_points={},slip=0,work=0) for _ in range(len(ends)+1)]
        for r in rows:
            band=sum(r["from"]>=t for t in ends)
            b=bins[band]
            w=r["to"]-r["from"]
            if r["accepted"]:
                prev_width=w
                b["accepted"]+=1
                continue
            b["rejected"]+=1
            b["grew_and_rejected"]+=int(prev_width is not None and w>prev_width)
            a=r["audit"] or {}
            if a.get("event"):
                key=str(a["event"][1])
                b["switch_points"][key]=b["switch_points"].get(key,0)+1
            b["slip"]+=int(bool(a) and a["max_slip"]>a["slip_budget"])
            b["work"]+=int(bool(a) and a["max_work"]>a["work_budget"])
        assert sum(b["accepted"]+b["rejected"] for b in bins)==len(rows)
        locations[label]=bins
    save(out/"cost-by-load.json",dict(parent_load_boundaries=ends,locations=locations))
    save(out/"comparison.json",dict(common=common,budget=budgets,unchanged_accepted_prefix=len(prior),
        switches=transitions,post_entrance_regrowth=[dict(until=t,groups=sum(g["accepted"]<=t for g in groups),
        rejections=sum(len(g["rejected"]) for g in groups if g["accepted"]<=t)) for t in (*FIXED,1.)]))


if __name__=="__main__":
    {"freeze":freeze,"score":score}[sys.argv[1]](Path(sys.argv[2]))
