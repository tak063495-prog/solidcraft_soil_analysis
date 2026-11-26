"""Freeze reload windows from archived observations, before solving."""
import gzip
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
P7 = HERE.parent / "2026-11-26-remaining-p7/c01"
P8 = HERE.parent / "2026-11-26-remaining-p8/c01"


def status(n):
    folder = P7 / f"original{n}"
    meta = json.loads((folder / "manifest.json").read_text())["fields"]["contact_status"]
    with zipfile.ZipFile(P7 / f"first-four-n{n}.zip") as pre, zipfile.ZipFile(folder / "fields-tail.zip") as tail:
        data = pre.read("contact_status.bin") + tail.read("contact_status.bin")
    assert hashlib.sha256(data).hexdigest() == meta["sha256"]
    values = np.frombuffer(data, dtype=meta["dtype"]).reshape(-1, 32)
    return values[5*256-1:6*256], meta["sha256"]


def events(n):
    values, digest = status(n)
    assert len(values) == 257
    events = {}
    for j, i in zip(*np.nonzero(values[:-1] != values[1:])):
        key = int(i), int(values[j,i]), int(values[j+1,i])
        events.setdefault(key, []).append(int(j))
    return events, digest


def freeze():
    coarse, a = events(256)
    medium, b = events(512)
    pairs, missing = [], []
    for key in sorted(coarse.keys() | medium.keys()):
        x, y = coarse.get(key, []), medium.get(key, [])
        if len(x) != len(y):
            missing.append(dict(key=key, coarse=x, medium=y))
            continue
        for j,k in zip(x,y,strict=True):
            pairs.append(dict(key=key, coarse_cell=j, medium_cell=k,
                              lower=min(j,k)/256, upper=(max(j,k)+1)/256))
    cells = sorted({j for p in pairs for j in range(int(p["lower"]*256),int(p["upper"]*256))})
    table = (P8 / "predictor.json").read_bytes()
    result = dict(status_hashes=[a,b], p8_table_sha256=hashlib.sha256(table).hexdigest(),
                  pairs=pairs, unmatched=missing, cells=cells,
                  caps=[c if j in cells else None for j,c in enumerate(json.loads(table)["caps"])],
                  progress_rule="Parent in half-open union of matched cells: apply P8 cap with ordinary native grid landing. Outside: ordinary proposal, including immediately after leaving. Native fine grid hits every cell boundary. Match by point, directed transition and chronological ordinal only for equal counts.",
                  comparison="Report all three paths at full budget and greatest mutually accepted native 1/1024 grid coordinate.")
    return result


if __name__ == "__main__":
    print(json.dumps(freeze(), indent=2))

