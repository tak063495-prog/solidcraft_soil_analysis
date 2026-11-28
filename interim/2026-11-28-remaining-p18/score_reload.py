"""Join pre-trial point forecasts to the independent, original trial verdict."""
from collections import Counter
import csv

import numpy as np
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
P17 = HERE.with_name("2026-11-28-remaining-p17")
sys.path.insert(0,str(P17))
from reload_score17 import FIXED, bins, counts, read


def check_forecast(note,law,base,minimum):
    start,previous,target = note["p18_parent"],note["previous"],note["nominal"]
    width = target-start
    predicted = []
    roots = []
    for i,((status,jump,normal,traction),old) in enumerate(zip(note["points"],note["previous_jumps"],strict=True)):
        jump,n,tau = map(np.array,(jump,normal,traction))
        delta = (jump-old)*(width/(start-previous))
        gap,slope = n@jump,n@delta
        def add(fraction,kind,guard):
            offset = width*fraction
            roots.append(dict(point=i,kind=kind,progress=start+offset,from_start=offset,from_end=width-offset,guard=guard))
            if guard < offset and guard < width-offset:
                predicted.append((i,kind,start+offset))
        if (gap>0)!=(gap+slope>0):
            guard = minimum if status=="slip" and slope>0 else base/64
            add(-gap/slope,"gap",guard)
        if status=="stick" and gap<=0:
            shear = tau-n*(n@tau)
            change = law["tangential_stiffness_kpa_m"]*(delta-n*slope)
            cap = law["friction_coefficient"]*law["normal_stiffness_kpa_m"]
            def margin(f):
                return np.linalg.norm(shear+change*f)+cap*(gap+slope*f)
            hi = -gap/slope if gap+slope>0 else 1.
            if margin(0)<0<margin(hi):
                lo = 0.
                for _ in range(48):
                    mid = (lo+hi)/2
                    if margin(mid)>0:
                        hi = mid
                    else:
                        lo = mid
                add((lo+hi)/2,"friction",base/64)
    for (i,kind,t),(j,method,to) in zip(predicted,note["forecasts"],strict=True):
        assert (i,kind)==(j,method)
        np.testing.assert_allclose(t,to,rtol=0,atol=4*np.spacing(target))
    expected = min([target]+[value for _,_,value in predicted])
    np.testing.assert_allclose(note["proposal"],expected,rtol=0,atol=4*np.spacing(target))
    return roots


def run(out):
    previous = json.loads((P17/"c01/reload/comparison.json").read_text())
    comparison = {k:previous[k] for k in ("common","budget","retrospective_comparator","prior_accepted_landing_table_used")}
    case=json.loads((HERE.with_name("2026-11-26-remaining-p7")/"c01/original1024/input.json").read_text())
    joint,=case["joint_system"]["joints"]
    base=next(s["initial_increment"] for s in case["stages"] if s["name"]=="reload")
    point_table = []
    for mode,folder in (("shadow",out),("control",out/"control")):
        trace = read(folder/"trace.jsonl.gz")
        if mode == "shadow":
            assert [r for r in trace if "p18_parent" not in r] == read(P17/"c01/reload/trace.jsonl.gz")
        rows = [r for r in trace if r.get("stage")=="reload" and "accepted" in r]
        notes = {r["attempt"]:r for r in trace if "p18_parent" in r}
        assert not any(r.get("stage")==5 and "charge" in r for r in trace)
        assert len(rows)<=1787 and (len(rows)==1787 or rows[-1]["to"]==1)
        last_parent = None
        states = {}
        stats = Counter()
        missed = []
        by_point = [Counter() for _ in range(32)]
        for i,r in enumerate(rows,1):
            start = r["from"]
            if .625<=start<.6875:
                n = notes.pop(i)
                assert n["p18_parent"]==start and n["previous"]==last_parent
                assert r["to"]==(n["nominal"] if mode=="shadow" else n["proposal"])
                roots = check_forecast(n,joint["material"],base,case["settings"]["minimum_increment"])
                if last_parent in states:
                    np.testing.assert_array_equal(n["previous_jumps"],[v[1] for v in states[last_parent]])
                if start in states:
                    assert states[start]==n["points"]
                states[start]=n["points"]
                forecasts = n["forecasts"]
                event = r["audit"]["event"]
                nearest = min(forecasts,key=lambda x:x[2],default=None)
                ids = {f[0] for f in forecasts}
                observed = None if event is None else event[1]
                stats["proposals"]+=1
                stats["forecast"]+=bool(forecasts)
                stats["event"]+=event is not None
                stats["forecast_without_event"]+=bool(forecasts) and event is None
                stats["same_point_in_set"]+=observed in ids
                stats["same_first_point"]+=nearest is not None and nearest[0]==observed
                if observed is not None and observed not in ids:
                    missed.append(dict(attempt=i,point=observed,parent=start,to=r["to"],event=event[0],
                                       roots=[x for x in roots if x["point"]==observed]))
                for j,point in enumerate(by_point):
                    point["forecast"] += j in ids
                    point["actual"] += j == observed
                    point["hit"] += j == observed and j in ids
                    point["false_alarm"] += j != observed and j in ids
                    point["miss"] += j == observed and j not in ids
                if event is not None or forecasts:
                    point_table.append(dict(mode=mode,attempt=i,parent=start,previous=n["previous"],
                        nominal=n["nominal"],forecast_points=json.dumps(forecasts,separators=(",",":")),
                        tried=r["to"],accepted=r["accepted"],event_point=observed,
                        event_at=None if event is None else event[0],same_point=observed in ids,
                        slip=r["audit"]["max_slip"]>r["audit"]["slip_budget"]))
            if r["accepted"]:
                a=r["audit"]
                assert r["force"]<=1e-10 and a["event"] is None
                assert a["max_slip"]<=a["slip_budget"] and a["max_work"]<=a["work_budget"]
                last_parent=start
            else:
                assert r["parent_unchanged"]
        assert not notes
        has_second = any(r["accepted"] and r["to"]==FIXED[-1] for r in rows)
        comparison[mode]=dict(point_forecast=stats,by_point=by_point,missed_events=missed,
                              bands=bins(rows)["width_1_over_16"] if has_second else None)
        if mode=="control":
            for t in FIXED:
                i=next((i for i,r in enumerate(rows,1) if r["accepted"] and r["to"]==t),None)
                comparison["common"][str(t)]["p18"]=counts(rows[:i]) if i else None
            comparison["budget"]["p18"]=counts(rows)
    with (out/"points.csv").open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(point_table[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(point_table)
    (out/"comparison.json").write_text(json.dumps(comparison,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
