"""Compare accepted switch-controlled full histories on the original binary coordinates."""
import argparse
from contextlib import ExitStack
import hashlib
import io
import json
from pathlib import Path
import sys

import numpy as np

from sealed import HERE, archive

sys.path.insert(0, str(HERE.parent / "2026-11-15-wet-joint-convergence"))
import study
sys.path.insert(0, str(HERE.parent / "2026-12-19-c04-transient-water-cause-diagnosis"))
from c04cause_metrics import cosine, metric

FIELDS = ("plastic_slip", "plastic_dissipation", "elastic_energy")


def series(read):
    manifest = json.loads(read("manifest.json"))
    rows = {}
    for chunk in manifest["chunks"]:
        raw = read(chunk["file"])
        if hashlib.sha256(raw).hexdigest() != chunk["sha256"]:
            raise ValueError("changed accepted fields")
        with np.load(io.BytesIO(raw), allow_pickle=False) as data:
            decoded = {k: study.codec.decode(data[k], manifest["encoding"])
                       if manifest.get("encoding") else data[k] for k in FIELDS}
            for i, row in enumerate(chunk["rows"]):
                key = row["stage"], row["progress"]
                if key in rows:
                    raise ValueError("duplicate accepted coordinate")
                rows[key] = {k: decoded[k][i].copy() for k in FIELDS}
    if len(rows) != manifest["count"]:
        raise ValueError("incomplete accepted fields")
    return rows


def compare(solutions, keys):
    if any(k not in s for s in solutions for k in keys):
        raise ValueError("missing accepted coordinate")
    result = {}
    for field in FIELDS:
        coarse, fine = [np.array([s[k][field] for k in keys]) for s in solutions]
        m = metric(coarse, fine)
        m["maximum_at"] = [*keys[m["maximum_index"][0]], *m["maximum_index"][1:]]
        result[field] = m
    return result


def measure(native):
    histories = {}
    with ExitStack() as stack:
        low, high = [stack.enter_context(archive(HERE.parent / name, "archive-index.json"))
                     for name in ("2026-12-11-c01-global-step-control", "2026-12-13-c01-cause-diagnosis")]
        for n in (128, 256, 512):
            histories["native", n] = series(lambda f: (native / f"g{n}/native" / f).read_bytes())
            z, prefix = (high, "experiments/switch-n512-q4") if n == 512 else (low, f"paths/g{n}")
            histories["independent", n] = series(lambda f: z.read(f"{prefix}/independent/{f}"))
    comparisons = {}
    direction = {}
    def keys(n):
        return [(stage, i/n) for stage in study.STAGES for i in range(1, n+1)]
    for n in (128, 256, 512):
        comparisons[f"native-independent/{n}"] = compare(
            [histories[side,n] for side in ("native", "independent")], keys(n))
    for side in ("native", "independent"):
        for n in (128, 256):
            comparisons[f"{side}/{n}-{2*n}"] = compare([histories[side,n], histories[side,2*n]], keys(n))
        direction[side] = {}
        for field in FIELDS:
            a,b,c = [np.array([histories[side,n][k][field] for k in keys(128)]) for n in (128,256,512)]
            n0, n1 = (np.linalg.norm(x) for x in (a-b,b-c))
            direction[side][field] = dict(coarse_norm=n0, fine_norm=n1, cosine=cosine(a-b,b-c),
                                          observed_order=float(np.log2(n0/n1)) if min(n0,n1)>0 else None)
    return dict(mesh_tet10=12, contact_gp=32, switch_bisections=4, comparisons=comparisons,
                counts={f"{s}/{n}":len(v) for (s,n),v in histories.items()}, directions=direction,
                extrapolated=None, gci_fine=None, continuum_error=None, passed=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("native", type=Path)
    p.add_argument("output", type=Path)
    a=p.parse_args()
    with a.output.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(measure(a.native), f, indent=2, allow_nan=False)
        f.write("\n")
