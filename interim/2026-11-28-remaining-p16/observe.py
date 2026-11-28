"""Sample measured recontact endpoints shared by all six archived histories."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P14 = HERE.with_name("2026-11-27-remaining-p14")
P15 = HERE.with_name("2026-11-28-remaining-p15")


def run(work, out, name, n):
    with gzip.open(P15/"c01/alternate/common-middle.jsonl.gz", "rt") as f:
        grid = {json.loads(s)["t"] for s in f}
    assert len(grid) == 382
    prior = (P15 if name == "outer" else P14)/"c01/alternate"/name/str(n)
    case = prior/"input.json"
    prefix = f"{name}-{n}"
    env = dict(os.environ, P16_FIELDS="1", P7_SWITCH_HALF="1", P11_OBSERVE="1",
               P10_KNOTS=str(P14/"c01/alternate"/name/"plan.json"))
    if name == "outer":
        env["P15_OUTER"] = "1"
    rows = []
    cmd = [work/"probe/target/release/examples/c01_methods","examples/joint-water-column-tet10-mesh.json",
           "examples/ground-joint-water-mc.json",case,work/(prefix+".json")]
    with (work/(prefix+".trace")).open("wb") as log, subprocess.Popen(
            cmd, env=env, stdout=subprocess.PIPE, stderr=log) as proc:
        for line in proc.stdout:
            r = json.loads(line)
            if r["t"] in grid:
                rows.append(r)
        if proc.wait():
            raise RuntimeError(proc.returncode)
    assert [r["t"] for r in rows] == sorted(grid)
    actual = json.loads((work/(prefix+".json")).read_text())
    archived = json.loads((prior/"result.json").read_text())
    assert actual["stages"] == archived["stages"]
    np.savez_compressed(out/(prefix+".npz"), t=sorted(grid), **{
        key:np.array([r[key] for r in rows]) for key in ("u","load")},
        **{key:np.array([[q[i] for q in r["points"]] for r in rows]) for i,key in
           enumerate(("traction","jump","plastic_slip","accumulated"))})


if __name__ == "__main__":
    run(Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3],int(sys.argv[4]))
