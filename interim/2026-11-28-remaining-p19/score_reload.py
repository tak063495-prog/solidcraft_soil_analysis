"""Match the band-9 forecasts to the solved trials and the two fixed loads."""
from collections import Counter
import csv
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P18 = HERE.with_name("2026-11-28-remaining-p18")
spec = importlib.util.spec_from_file_location("score18",P18/"score_reload.py")
p18 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p18)


def run(out):
    previous = json.loads((P18/"c01/reload/comparison.json").read_text())
    result = {k:previous[k] for k in ("common","budget","retrospective_comparator","prior_accepted_landing_table_used")}
    case = json.loads((HERE.with_name("2026-11-26-remaining-p7")/"c01/original1024/input.json").read_text())
    base = next(s["initial_increment"] for s in case["stages"] if s["name"]=="reload")
    law = case["joint_system"]["joints"][0]["material"]
    points = []
    for mode,d in (("shadow",out),("control",out/"control")):
        trace = p18.read(d/"trace.jsonl.gz")
        if mode=="shadow":
            legacy = p18.read(P18/"c01/reload/control/trace.jsonl.gz")
            assert [r for r in trace if not ("p18_parent" in r and r["p18_parent"]<.625)] == legacy
        notes = {r["attempt"]:r for r in trace if "p18_parent" in r}
        rows = [r for r in trace if r.get("stage")=="reload" and "accepted" in r]
        assert not any(r.get("stage")==5 and "charge" in r for r in trace)
        assert len(rows)<=1787
        by_point = [Counter() for _ in range(32)]
        stats = Counter()
        reasons = Counter()
        states = {}
        previous_parent = None
        for i,r in enumerate(rows,1):
            start = r["from"]
            audit = r["audit"]
            if .5625 <= start < .6875:
                n = notes[i]
                assert n["previous"]==previous_parent and n["p18_parent"]==start
                if start in states:
                    assert states[start] == n["points"]
                if previous_parent in states:
                    np.testing.assert_array_equal(n["previous_jumps"],[v[1] for v in states[previous_parent]])
                states[start] = n["points"]
                eligible = [f for f in n["forecasts"] if f[0]!=19]
                proposed = min([n["nominal"]]+[f[2] for f in eligible])
                assert proposed==n["proposal"]
                shortened = start>=.625 or mode=="control"
                assert r["to"]==(proposed if shortened else n["nominal"])
                if shortened and proposed<n["nominal"]:
                    assert all(f[0]!=19 for f in eligible if f[2]==proposed)
                # The p18 arithmetic checker validates each point and both exclusion guards.
                observed_note = dict(n,proposal=min([n["nominal"]]+[f[2] for f in n["forecasts"]]))
                roots = p18.check_forecast(observed_note,law,base,case["settings"]["minimum_increment"])
                if start<.625:
                    forecasts = n["forecasts"]
                    event = audit["event"]
                    actual = event[1] if event else None
                    ids = {f[0] for f in forecasts}
                    first = min(forecasts,key=lambda f:f[2],default=None)
                    slip = audit["max_slip"]>audit["slip_budget"]
                    work = audit["max_work"]>audit["work_budget"]
                    label = "+".join(k for k,v in (("event",event),("slip",slip),("work",work)) if v)
                    reasons["accepted" if r["accepted"] else label or r["reason"]] += 1
                    stats["proposals"] += 1
                    stats["with_forecast"] += bool(forecasts)
                    stats["with_event"] += bool(event)
                    stats["same_first_point"] += first is not None and first[0]==actual
                    stats["same_point_in_set"] += actual in ids
                    stats["forecast_without_event"] += bool(forecasts) and not event
                    stats["shortened"] += r["to"]<n["nominal"]
                    stats["shortened_accepted"] += r["to"]<n["nominal"] and r["accepted"]
                    for j,c in enumerate(by_point):
                        c["forecast"] += j in ids
                        c["actual"] += actual==j
                        c["hit"] += actual==j and j in ids
                        c["miss"] += actual==j and j not in ids
                    if forecasts or event or not r["accepted"]:
                        points.append(dict(mode=mode,attempt=i,parent=start,nominal=n["nominal"],tried=r["to"],
                            forecasts=json.dumps(forecasts),all_roots=json.dumps(roots),event_at=event[0] if event else None,
                            event_point=actual,slip=slip,accepted=r["accepted"]))
            if r["accepted"]:
                assert audit["event"] is None and r["force"]<=1e-10
                assert audit["max_slip"]<=audit["slip_budget"] and audit["max_work"]<=audit["work_budget"]
                previous_parent=start
            else:
                assert r["parent_unchanged"]
        summary = json.loads((d/"result.json").read_text())
        prefix = json.loads((out/"prefix-identity.json").read_text())["prefix_charges"]
        stages = {}
        for j,s in enumerate(summary["stages"]):
            extra = sum(r.get("stage")==j and "charge" in r for r in trace)
            if s["stage"] in prefix:
                extra = prefix[s["stage"]]-s["attempts"]
            assert extra>=0
            charged = s["attempts"]+extra
            assert charged<=1787
            stages[s["stage"]] = dict(outer=s["attempts"],extra=extra,charged=charged,progress=s["last_accepted"])
        result[mode] = dict(band9=stats,reasons=reasons,points=by_point,stages=stages,
                            complete=summary["complete"],stop=summary["stop"])
        if mode=="control":
            for t in p18.FIXED:
                i = next((i for i,r in enumerate(rows,1) if r["accepted"] and r["to"]==t),None)
                result["common"][str(t)]["p19"] = p18.counts(rows[:i]) if i else None
            result["budget"]["p19"] = p18.counts(rows)
    with (out/"band9.csv").open("w") as f:
        writer = csv.DictWriter(f,fieldnames=list(points[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(points)
    (out/"comparison.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
