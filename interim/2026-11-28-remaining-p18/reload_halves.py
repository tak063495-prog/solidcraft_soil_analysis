"""Replay the existing six histories and keep only the accepted reload half-intervals."""
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


def replay(work, out, route, n):
    archive = (P15 if route == "outer" else P14) / "c01/alternate" / route / n
    env = dict(os.environ, P18_HALVES="1", P7_SWITCH_HALF="1", P11_OBSERVE="1",
               P10_KNOTS=str(P14/"c01/alternate"/route/"plan.json"))
    if route == "outer":
        env["P15_OUTER"] = "1"
    name = f"{route}-{n}"
    result = work / (name + ".json")
    trace = work / (name + ".trace")
    binary = work / "probe/target/release/examples/c01_methods"
    cmd = [binary, "examples/joint-water-column-tet10-mesh.json",
           "examples/ground-joint-water-mc.json", archive / "input.json", result]
    candidates = {}
    with trace.open("wb") as err, subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, env=env) as proc:
        for line in proc.stdout:
            row = json.loads(line)
            start, _, end = row["span"]
            assert (start, end) not in candidates
            candidates[start, end] = row
        if proc.wait():
            raise RuntimeError(proc.returncode)
    assert json.loads(result.read_text())["stages"] == json.loads((archive/"result.json").read_text())["stages"]
    with trace.open("rb") as current, gzip.open(archive/"trace.jsonl.gz", "rb") as saved:
        digest = hashlib.file_digest(current, "sha256").hexdigest()
        assert digest == hashlib.file_digest(saved, "sha256").hexdigest()
    with trace.open() as stream:
        accepted = [r for line in stream if (r := json.loads(line)).get("stage") == "reload" and r["accepted"]]
    t, states, spans = [], [], []
    for r in accepted:
        row = candidates[r["from"], r["to"]]
        start, mid, end = row["span"]
        if not t:
            assert start == 0
            t.append(start)
            states.append(row["initial"])
            normals = row["normals"]
        assert start == t[-1] and start < mid < end
        t.extend([mid, end])
        states.extend([row["half"], row["end"]])
        spans.append(row["span"])
    assert t[-1] == 1
    np.savez_compressed(work/(name+".npz"), t=t, spans=spans, normals=normals,
                        **{k: [r[k] for r in states] for k in ("u", "load", "traction", "jump")})
    (out/(name+".json")).write_text(json.dumps(dict(trace_sha256=digest, accepted=len(spans),
        endpoints=len(t), stages=json.loads(result.read_text())["stages"]), indent=2)+"\n")


if __name__ == "__main__":
    replay(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), *sys.argv[3:])
