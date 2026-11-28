"""Collect native half and fine states, filtering against the final accepted ledger."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P14 = HERE.with_name("2026-11-27-remaining-p14")
P15 = HERE.with_name("2026-11-28-remaining-p15")


def run(work, out, route, n):
    prior = (P15 if route=="outer" else P14)/"c01/alternate"/route/str(n)
    name = f"{route}-{n}"
    env = dict(os.environ,P17_RETURN="1",P7_SWITCH_HALF="1",P11_OBSERVE="1",
               P10_KNOTS=str(P14/"c01/alternate"/route/"plan.json"))
    if route=="outer":
        env["P15_OUTER"]="1"
    result = work/(name+".json")
    cmd = [work/"probe/target/release/examples/c01_methods","examples/joint-water-column-tet10-mesh.json",
           "examples/ground-joint-water-mc.json",prior/"input.json",result]
    candidates = {}
    trace = work/(name+".trace")
    with trace.open("wb") as log, subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log,env=env) as proc:
        for line in proc.stdout:
            row = json.loads(line)
            start,mid,end = row["span"]
            assert (start,end) not in candidates
            candidates[start,end] = row
        if proc.wait():
            raise RuntimeError(proc.returncode)
    assert json.loads(result.read_text())["stages"] == json.loads((prior/"result.json").read_text())["stages"]
    with gzip.open(prior/"trace.jsonl.gz","rb") as archived,trace.open("rb") as current:
        assert hashlib.file_digest(current,"sha256").digest()==hashlib.file_digest(archived,"sha256").digest()
    with trace.open() as stream:
        accepted = [r for line in stream if (r:=json.loads(line)).get("stage")=="return" and r["accepted"]]
    series = []
    times = []
    spans = []
    for r in accepted:
        s = candidates[r["from"],r["to"]]
        start,mid,end = s["span"]
        if not times:
            assert start==0
            normals = s["normals"]
            times.append(start)
            series.append(s["initial"])
        assert times[-1]==start
        spans.append([start,mid,end])
        times.extend([mid,end])
        series.extend([s["half"],s["end"]])
    assert times[-1]==1 and np.all(np.diff(times)>0)
    np.savez_compressed(out/(name+".npz"),t=times,spans=spans,normals=normals,**{
        key:np.array([s[key] for s in series]) for key in ("u","load","traction","jump")})


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve(),sys.argv[3],int(sys.argv[4]))
