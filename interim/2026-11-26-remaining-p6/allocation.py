"""Same-parent diagnostic comparisons on the common prefix of two actual routes."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P5 = HERE.parent/"2026-11-26-remaining-p5"
sys.path.insert(0,str(P5))
import c01_trace as old


def read(path):
    with gzip.open(path,"rt") as stream:
        return list(map(json.loads,stream))


def choose(rows):
    """Accepted switches for every point, plus largest accepted/rejected errors nearby."""
    picked = set()
    for idx, r in enumerate(rows):
        statuses = [r[k]["status"] for k in ("parent","after_m","after_f")]
        if r["accepted"] and 234 < 256*r["end"] <= 239 and any(a!=b for a,b in zip(statuses,statuses[1:])):
            picked.add(idx)
    for lo, hi in ((234,235),(235,236),(236,237),(237,239),(239,240)):
        for accepted in (True,False):
            region = [i for i,r in enumerate(rows) if lo < r["end"]*256 <= hi and r["accepted"]==accepted]
            if region:
                picked.add(max(region,key=lambda i:rows[i]["post_max_slip"]))
    return sorted(picked)


def run(out, binary, n):
    out.mkdir(parents=True)
    rows = read(P5/f"c01/n{n}/trials.jsonl.gz")
    chosen = choose(rows)
    audit = json.loads((P5/f"c01/n{n}/audit.json").read_text())
    spent = audit["attempts"]+audit["additional_material_and_global_trials"]
    assert spent + 12*len(chosen) <= 2065
    windows = [[rows[i]["start"],rows[i]["end"]] for i in chosen]
    assert len({tuple(w) for w in windows})==len(windows)
    case = json.loads((P5/f"c01/n{n}/prefix.json").read_bytes())
    for src in (old.P4/f"original{n}-input.json",old.P3/f"case-{n}.json"):
        full = json.loads(src.read_bytes())
        full["stages"] = full["stages"][:3]
        assert full == case
    cmd = [binary,old.ROOT/"examples/joint-water-column-tet10-mesh.json",
           old.ROOT/"examples/ground-joint-water-mc.json",P5/f"c01/n{n}/prefix.json",out/"result.json"]
    actual = []
    matched = 0
    with subprocess.Popen(cmd,stderr=subprocess.PIPE,env=dict(os.environ,P6_WINDOWS=json.dumps(windows))) as proc, \
            gzip.open(out/"candidates.jsonl.gz","wt") as sink, gzip.open(out/"parents.jsonl.gz","wt") as parents:
        for line in proc.stderr:
            try:
                r=json.loads(line)
            except ValueError:
                raise RuntimeError(line.decode(errors="replace"))
            if not r.get("allocation"):
                assert r==rows[matched]
                matched+=1
                continue
            idx=chosen[len(actual)]
            ref=rows[idx]
            assert [r["start"],r["end"]]==windows[len(actual)] and ref["accepted"]==r["accepted"]
            for a,b in (("parent","parent"),("coarse","after_c"),("half","after_m"),("fine","after_f")):
                assert all(r[a][key]==ref[b][key] for key in ("status","slip"))
            parents.write(json.dumps(dict(index=idx,parent=r.pop("saved_parent")))+"\n")
            r["index"]=idx
            sink.write(json.dumps(r)+"\n")
            actual.append(r["extra_solves"])
        if proc.wait():
            raise RuntimeError("diagnostic replay failed")
    assert len(actual)==len(chosen) and matched==len(rows)
    baseline=json.loads((P5/f"c01/n{n}/result.json").read_bytes())
    route=json.loads((out/"result.json").read_bytes())
    assert route["complete"] and route["stages"]==baseline["stages"]
    assert spent+sum(actual)<=2065
    (out/"audit.json").write_text(json.dumps(dict(n=n, original_and_alternate_prefix_equal=True,
        matched_p5_trials=matched,selected=chosen,base_attempts=audit["attempts"],
        p5_extra=spent-audit["attempts"],extra_solve_calls=sum(actual),charged=spent+sum(actual),
        cap=2065),indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]),Path(sys.argv[2]).resolve(),int(sys.argv[3]))
