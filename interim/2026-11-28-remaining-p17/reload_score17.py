"""Match pre-trial proposals to their latest accepted audit and untouched stage verdict."""
from collections import Counter
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
P16=HERE.with_name("2026-11-28-remaining-p16")
sys.path.insert(0,str(P16))
from reload_cost import FIXED,bins,counts,read


def run(out):
    trace=read(out/"trace.jsonl.gz")
    rows=[r for r in trace if r.get("stage")=="reload" and "accepted" in r]
    assert not any(r.get("stage")==5 and "charge" in r for r in trace)
    assert len(rows)<=1787 and (rows[-1]["to"]==1 or len(rows)==1787)
    comp=json.loads((P16/"c01/reload/comparison.json").read_text())
    comp={k:comp[k] for k in ("common","budget","retrospective_comparator","prior_accepted_landing_table_used")}
    rules=[("p16_parent","last_slip_ratio",.625,.6875,"max_slip","slip_budget"),
           ("p17_parent","last_work_ratio",.8125,.875,"max_work","work_budget")]
    for tag,field,lo,hi,err,budget in rules:
        notes={r["attempt"]:r for r in trace if tag in r}
        proposals=len(notes)
        hits=Counter()
        last=None
        for i,r in enumerate(rows,1):
            p=r["from"]
            if lo<=p<hi:
                note=notes.pop(i)
                assert note[tag]==p and last==(note["last_width"],note[field])
                h=note["nominal"]-p
                value=last[1]*(h/last[0])**2
                assert abs(value-note["predicted_ratio"])<=2e-15*max(1,value)
                assert note["to"]==(p+h/2 if value>1 else note["nominal"])
                # The inherited downstream event predictor may further shorten this proposal.
                assert r["to"]<=note["to"]
                if value>1:
                    hits["accepted" if r["accepted"] else "rejected"]+=1
            if r["accepted"]:
                a=r["audit"]
                assert r["force"]<=1e-10 and a["event"] is None
                assert a["max_work"]<=a["work_budget"] and a["max_slip"]<=a["slip_budget"]
                last=r["to"]-p,a[err]/a[budget]
            else:
                assert r["parent_unchanged"]
        assert not notes
        comp[tag]=dict(proposals=proposals,shortened=hits)
    for t in FIXED:
        i=next((i for i,r in enumerate(rows,1) if r["accepted"] and r["to"]==t),None)
        comp["common"][str(t)]["p17"]=counts(rows[:i]) if i else None
    comp["budget"]["p17"]=counts(rows)
    comp["band_counts_to_second_load"]=bins(rows)["width_1_over_16"]
    comp["reload_finished"]=comp["budget"]["p17"]["progress"]==1
    (out/"comparison.json").write_text(json.dumps(comp,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
