import gzip
import json
from pathlib import Path
import sys


def load_rows(path):
    with gzip.open(path,"rt") as f:
        return [json.loads(line) for line in f if line.startswith("{")]


def run(folder):
    p8=Path(__file__).resolve().parent.parent/"2026-11-26-remaining-p8/c01"
    series={}
    result={}
    for label in ("baseline","predicted","windows"):
        rows=[r for r in load_rows(folder/f"{label}-trace.jsonl.gz") if r.get("stage")=="reload"]
        assert all(r["parent_unchanged"] for r in rows if not r["accepted"])
        if label!="windows":
            prior=[r for r in load_rows(p8/f"{label}-trace.jsonl.gz") if r.get("stage")=="reload"]
            assert [(r["from"],r["to"],r["accepted"],r["reason"]) for r in rows]==[
                   (r["from"],r["to"],r["accepted"],r["reason"]) for r in prior]
        series[label]=rows
        result[label]=dict(attempts=len(rows),accepted=sum(r["accepted"] for r in rows),
                           rejected=sum(not r["accepted"] for r in rows),
                           end=max(r["to"] for r in rows if r["accepted"]))
        assert len(rows)<=1787
    windows=json.loads((Path(__file__).parent/"windows.json").read_text())["caps"]
    def target(parent,width):
        import math,sys
        base=1/1024
        boundary=(math.floor(parent/base)+1)*base
        if boundary<=parent+16*sys.float_info.epsilon: boundary+=base
        return min(1.,boundary,parent+width)
    audits={}
    p8caps=json.loads((p8/"predictor.json").read_text())["caps"]
    for label,rows in series.items():
        inc=1/1024
        shortened=outside=0
        for r in rows:
            parent=r["from"]
            proposal=target(parent,inc)
            cap=windows[int(parent*256)] if label=="windows" else None
            if label=="predicted":
                cap=p8caps[int(parent*256)]
            expected=min(proposal,target(parent,cap)) if cap is not None and proposal-parent>cap else proposal
            assert expected==r["to"],(label,parent,expected,r["to"])
            outside+=int(cap is None)
            shortened+=int(expected<proposal)
            audit=r["audit"]
            inc=audit["proposed_increment"] if audit and audit["proposed_increment"] is not None else (r["to"]-parent)/2
        audits[label]=dict(shortened=shortened,ordinary_region=outside,verified=len(rows))
    shared=set.intersection(*[{r["to"] for r in rows if r["accepted"] and r["to"]*1024==int(r["to"]*1024)}
                              for rows in series.values()])
    def matched(t):
        out={}
        for label,rows in series.items():
            last=next(i for i,r in enumerate(rows) if r["accepted"] and r["to"]==t)+1
            out[label]=dict(attempts=last,rejected=sum(not r["accepted"] for r in rows[:last]))
        return dict(t=t,paths=out)
    (folder/"comparison.json").write_text(json.dumps(dict(full_budget=result,proposal_audit=audits,matched=matched(max(shared)),
        grid_comparisons=[matched(t) for t in sorted(shared)]),indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
