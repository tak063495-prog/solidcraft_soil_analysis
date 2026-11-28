"""Count mutually exclusive rejection reasons and replay only the proposal arithmetic."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
P15 = HERE.with_name("2026-11-28-remaining-p15")
sys.path.insert(0,str(HERE.with_name("2026-11-27-remaining-p14")))
from reload_score import FIXED,counts,read


def bins(rows):
    cut = next(i for i,r in enumerate(rows,1) if r["accepted"] and r["to"]==FIXED[-1])
    bands = defaultdict(Counter)
    grid = defaultdict(Counter)
    for r in rows[:cut]:
        if r["accepted"]:
            label = "accepted"
        else:
            a = r["audit"] or {}
            reasons = []
            if a.get("event") is not None:
                reasons.append("event")
            if a and a["max_slip"]>a["slip_budget"]:
                reasons.append("slip")
            if a and a["max_work"]>a["work_budget"]:
                reasons.append("work")
            label = "+".join(reasons) or r["reason"]
            assert r["parent_unchanged"]
        for n,b in [(16,bands),(1024,grid)]:
            b[int(r["from"]*n)][label] += 1
    return dict(width_1_over_16=bands,width_1_over_1024=grid,common=counts(rows[:cut]))


def run(out):
    prior = read(P15/"c01/reload/trace.jsonl.gz")
    p15 = [r for r in prior if r.get("stage")=="reload" and "accepted" in r]
    a = bins(p15)
    assert (a["common"]["accepted"],a["common"]["rejected"])==(1627,99)
    (out/"p15-bands.json").write_text(json.dumps(a,indent=2)+"\n")
    if not (out/"trace.jsonl.gz").is_file():
        return
    trace = read(out/"trace.jsonl.gz")
    assert not any(x.get("stage")==5 and "charge" in x for x in trace)
    rows = [r for r in trace if r.get("stage")=="reload" and "accepted" in r]
    notes = {r["attempt"]:r for r in trace if "p16_parent" in r}
    last = None
    hits = Counter()
    for i,r in enumerate(rows,1):
        p = r["from"]
        audit = r["audit"]
        if .625<=p<.6875:
            n = notes.pop(i)
            assert last is not None and n["p16_parent"]==p
            assert last==(n["last_width"],n["last_slip_ratio"])
            h = n["nominal"]-p
            ratio = last[1]*(h/last[0])**2
            assert abs(n["predicted_ratio"]-ratio)<=2e-15*max(1,ratio)
            assert r["to"]==n["to"]==(p+h/2 if ratio>1 else n["nominal"])
            if ratio>1:
                hits["accepted" if r["accepted"] else "rejected"]+=1
        if r["accepted"]:
            assert r["force"]<=1e-10 and audit["event"] is None
            assert audit["max_work"]<=audit["work_budget"] and audit["max_slip"]<=audit["slip_budget"]
            last = (r["to"]-p,audit["max_slip"]/audit["slip_budget"])
        else:
            assert r["parent_unchanged"]
    assert not notes and len(rows)==1787
    comp = json.loads((P15/"c01/reload/comparison.json").read_text())
    comp = {k:comp[k] for k in ("common","budget","retrospective_comparator","prior_accepted_landing_table_used")}
    for t in FIXED:
        i = next((i for i,r in enumerate(rows,1) if r["accepted"] and r["to"]==t),None)
        comp["common"][str(t)]["p16"] = counts(rows[:i]) if i else None
    comp["budget"]["p16"] = counts(rows)
    comp.update(shortened=hits,p16_bands=bins(rows) if comp["common"][str(FIXED[-1])]["p16"] else None)
    (out/"comparison.json").write_text(json.dumps(comp,indent=2)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]))
