"""Freeze the reload controls from the committed p10 trace before new solves."""
import collections
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

def freeze(trace, plan):
    table = json.loads(plan.read_bytes())
    bands = table["windows"]
    events = [e for e in table["events"] if e["used"]]
    rows = []
    controls = set()
    increment = 1 / 1024
    attempts = 0
    with gzip.open(trace, "rt") as stream:
        for line in stream:
            if not line.startswith("{"):
                continue
            r = json.loads(line)
            if r.get("stage") != "reload":
                continue
            attempts += 1
            p = r["from"]
            knot = (math.floor(p*1024)+1)/1024
            if knot <= p + 16*sys.float_info.epsilon:
                knot += 1/1024
            normal = min(1., knot, p+increment)
            outside = not any(l <= p < h for l,h in bands)
            if outside:
                assert normal == r["to"]
                hits = [dict(span=e["span"], point=e["point"], transition=e["transition"],
                             complete=normal>=e["span"][1]) for e in events if p<e["span"][0]<normal]
                if hits:
                    rows.append(dict(attempt=attempts, parent=p, normal=normal, accepted=r["accepted"],
                                     reason=r["reason"], hits=hits))
                    if not r["accepted"]:
                        controls.update((tuple(h["span"]),tuple(h["transition"])) for h in hits if h["complete"])
            audit = r["audit"]
            increment = audit["proposed_increment"] if audit and audit["proposed_increment"] is not None else (r["to"]-p)/2
    counts = collections.defaultdict(collections.Counter)
    for row in rows:
        for key in {(" -> ".join(h["transition"]), h["complete"]) for h in row["hits"]}:
            counts[str(key)]["accepted" if row["accepted"] else "rejected"] += 1
    return dict(attempts=attempts, trace_sha256=hashlib.sha256(trace.read_bytes()).hexdigest(),
                plan_sha256=hashlib.sha256(plan.read_bytes()).hexdigest(), outside_entered=len(rows),
                outside_crossed=sum(any(h["complete"] for h in r["hits"]) for r in rows),
                rejected_crossed=sum(not r["accepted"] and any(h["complete"] for h in r["hits"]) for r in rows),
                by_direction=counts, rows=rows, bands=bands, controls=[dict(span=s, transition=t,
                target=(s[0]+s[1])/2 if t==("stick","slip") else s[0]) for s,t in sorted(controls)],
                rule="From outside all bands, clip a full forward crossing of a listed span/transition to its frozen target (stick-to-slip midpoint, slip-to-open left endpoint), taking the earliest target. Otherwise use ordinary control.")
if __name__ == "__main__":
    print(json.dumps(freeze(*map(Path,sys.argv[1:])),indent=2,allow_nan=False))
