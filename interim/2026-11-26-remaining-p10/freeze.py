"""Construct opt-in control data exclusively from prior accepted observations."""
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P9 = HERE.parent / "2026-11-26-remaining-p9"


def rows(path):
    with gzip.open(path, "rt") as stream:
        return [json.loads(line) for line in stream if line.startswith("{")]


def freeze_alt():
    h = np.load(P9 / "c01/alternate/accepted-halves.npz")
    clock = h["recontact_512_progress"]
    grid = json.loads((P9 / "c01/alternate/matched-plan.json").read_text())
    spans = []
    knots = set()
    for i, point in enumerate(grid["points"]):
        status = h["recontact_512_classes"][:, i]
        j, = np.flatnonzero(status[:-1] != status[1:])
        assert "open" in status[j].lower() and status[j+1]=="slip"
        l, r = clock[j:j+2]
        spans.append(dict(point=point, left=l, right=r, attempt_half=h["recontact_512_attempt_half"][j+1].tolist()))
        knots.update(np.linspace(l, r, 5).tolist())
        knots.update(t for t in grid["grids"][i] if l <= t <= r)
    return dict(spans=spans, knots=sorted(knots))


def freeze_reload(path):
    log = rows(path)
    actual = [r for r in log if r.get("stage") == "reload"]
    old = rows(P9 / "c01/reload/windows-trace.jsonl.gz")
    fields = ["from", "to", "accepted", "reason"]
    signature = lambda r: [r[k] for k in fields]
    assert list(map(signature, actual)) == [signature(r) for r in old if r.get("stage") == "reload"]
    candidates = {(r["span"][0], r["span"][-1]): r for r in log if "span" in r}
    bands = json.loads((P9 / "windows.json").read_text())["cells"]
    pairs = []
    windows = set()
    for index, r in enumerate(actual, 1):
        if not r["accepted"]:
            continue
        c = candidates[r["from"], r["to"]]
        for j in range(2):
            l, h = c["span"][j:j+2]
            for i, (a, b) in enumerate(zip(*c["labels"][j:j+2], strict=True)):
                if a != b:
                    used = int(l*256) in bands
                    pairs.append(dict(point=i, transition=[a,b], attempt=index, half=j+1, span=[l,h], used=used))
                    if used:
                        windows.add((l,h))
    assert windows
    return dict(trace_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), events=pairs, windows=sorted(windows))


if __name__ == "__main__":
    obj = freeze_alt() if len(sys.argv)==1 else freeze_reload(Path(sys.argv[1]))
    print(json.dumps(obj, indent=2, allow_nan=False))
