"""Score fixed-parent histories without sending a candidate back to the controller."""
import json
from pathlib import Path
import sys

import numpy as np

from allocation import P5, read


def region(t):
    if t <= 235:
        return "before"
    if t <= 236:
        return "switch_235_236"
    return "after"


def norm(a):
    return float(np.max(np.abs(a)))


def run(out):
    report = {}
    for n in (256,512,1024):
        folder = out/f"n{n}"
        rows = read(P5/f"c01/n{n}/trials.jsonl.gz")
        candidates = read(folder/"candidates.jsonl.gz")
        regs = {}
        for label in ("before","switch_235_236","after"):
            series = [r for r in rows if 234 < r["end"]*256 <= 240 and region(r["end"]*256)==label]
            reasons = {}
            changes = [set(),set()]
            max_ratio = 0.
            fixed_only = 0
            for r in series:
                s = r["post_max_slip"] > r["slip_budget"]
                d = r["post_max_dissipation"] > r["dissipation_budget"]
                e = r["event"] is not None
                assert r["accepted"] == (not(s or d or e))
                reason = "+".join(key for key,x in (("slip",s),("work",d),("event",e)) if x) or "accepted"
                reasons[reason] = reasons.get(reason,0)+1
                max_ratio = max(max_ratio,r["post_max_slip"]/r["slip_budget"])
                fixed_only += bool(s and not(d or e) and r["post_max_slip"] <= 1e-4/n**2)
                if r["accepted"]:
                    for i,(a,b,c) in enumerate(zip(*(r[k]["status"] for k in ("parent","after_m","after_f")))):
                        if (a,b)==("open","slip"): changes[0].add(i)
                        if (b,c)==("open","slip"): changes[1].add(i)
            regs[label] = dict(trials=len(series),reasons=reasons, fixed_only_vs_scaled=fixed_only,
                max_slip_budget_ratio=max_ratio,accepted_open_slip_by_half=[sorted(s) for s in changes],
                corrections_coarse_mid_end=[sum(r["selected_from"][j] is not None for r in series) for j in range(3)])
        scores = []
        for r in candidates:
            def diff(a,b,k):
                return norm(np.array(r[a][k])-r[b][k])
            row = dict(index=r["index"],start256=r["start"]*256,end256=r["end"]*256,
                       baseline_accepted=r["accepted"],event=r["event"],paths={})
            # Two local step-doubling defects on their own immutable half-parents.
            # One fraction beta of each original allowance goes to the first half.
            mid = {k:diff("early_mid","half",k) for k in ("slip","work")}
            end = {k:diff("late","fine",k) for k in ("slip","work")}
            bounds = {"slip":r["slip_budget"], "work":r["work_budget"]}
            lo = max(mid[k]/bounds[k] for k in bounds)
            hi = 1-max(end[k]/bounds[k] for k in bounds)
            row["half_defects"] = dict(mid=mid, end=end, min_mid_fraction=lo,
                max_mid_fraction=hi, equal_split_fits=lo<=.5<=hi,
                an_allocation_fits=lo<=hi)
            for name in ("fine","early","late"):
                slip = diff("coarse",name,"slip")
                work = diff("coarse",name,"work")
                row["paths"][name]=dict(coarse_slip_error=slip, coarse_work_error=work,
                    vs_two_half_slip=diff("fine",name,"slip"),vs_two_half_work=diff("fine",name,"work"),
                    within_both_local_bounds=slip<=r["slip_budget"] and work<=r["work_budget"])
            scores.append(row)
        stats={}
        for label in regs:
            rs=[r for r in scores if region(r["end256"])==label]
            parts=[r["half_defects"] for r in rs]
            stats[label]={"fixed_total_split":dict(equal_fits=sum(x["equal_split_fits"] for x in parts),
                feasible=sum(x["an_allocation_fits"] for x in parts),samples=len(parts),
                needs_more_mid=sum(x["an_allocation_fits"] and x["min_mid_fraction"]>.5 for x in parts),
                needs_more_end=sum(x["an_allocation_fits"] and x["max_mid_fraction"]<.5 for x in parts))}
            for name in ("early","late"):
                pairs=[(r["paths"]["fine"],r["paths"][name]) for r in rs]
                stats[label][name]=dict(samples=len(rs),
                    from_over_to_within=sum(not a["within_both_local_bounds"] and b["within_both_local_bounds"] for a,b in pairs),
                    from_within_to_over=sum(a["within_both_local_bounds"] and not b["within_both_local_bounds"] for a,b in pairs),
                    within_both=sum(b["within_both_local_bounds"] for a,b in pairs),
                    max_change_from_two_half_m=max(b["vs_two_half_slip"] for a,b in pairs))
        report[n]=dict(regions=regs,allocation_summary=stats,candidates=scores)
    (out/"summary.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
