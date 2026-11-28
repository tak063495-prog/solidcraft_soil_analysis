"""Pointwise classifications of all pre-trial forecasts, including false alarms."""
from collections import Counter
import csv
import json
from pathlib import Path
import sys


def run(out):
    info=json.loads((out/"comparison.json").read_text())
    trials=info["trials"]
    rejected=[]
    points=[Counter() for _ in range(32)]
    for r in trials:
        for i,counts in enumerate(points):
            predicted,actual=i in r["predicted"],i in r["actual"]
            key=("true_" if actual==predicted else "false_")+("positive" if predicted else "negative")
            counts[key]+=1
        if not r["accepted"]:
            p,a=r["predicted_argmax"],r["actual_argmax"]
            rejected.append(dict(attempt=r["attempt"],parent=r["parent"],to=r["to"],
                predicted_points=" ".join(map(str,r["predicted"])),actual_points=" ".join(map(str,r["actual"])),
                predicted_argmax=p,actual_argmax=a,predicted_max_ratio=r["predicted_ratio"][p],
                predicted_ratio_at_actual=r["predicted_ratio"][a],actual_max_ratio=r["actual_ratio"][a],
                exact_point_set=r["predicted"]==r["actual"]))
    with (out/"slip-rejections.csv").open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rejected[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rejected)
    data=dict(trials=len(trials),rejected=len(rejected),
        predicted_rejections=sum(bool(r["predicted"]) for r in trials),
        detected=sum(bool(r["predicted"]) and not r["accepted"] for r in trials),
        false_alarms=sum(bool(r["predicted"]) and r["accepted"] for r in trials),
        exact_rejected_sets=sum(r["exact_point_set"] for r in rejected),
        actual_peak_predicted=sum(r["predicted_ratio_at_actual"]>1 for r in rejected),
        peak_identity=sum(r["predicted_argmax"]==r["actual_argmax"] for r in rejected),
        by_point=[dict(point=i,**c) for i,c in enumerate(points)])
    (out/"slip-summary.json").write_text(json.dumps(data,indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
