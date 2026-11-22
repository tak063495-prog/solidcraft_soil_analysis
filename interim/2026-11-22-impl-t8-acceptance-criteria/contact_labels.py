"""Strict endpoint categories, last substep and accumulated work for the 13-stage cycle."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import numpy as np

from sealed import HERE, archive

SOURCE = HERE.parent / "2026-12-19-c04-transient-water-cause-diagnosis"
sys.path.insert(0, str(SOURCE))
import c04cause_metrics as metrics
import c04cause_observe as observe


def compare(root, mechanical, water):
    folders = [root / f"{kind}-m{mechanical}-w{water}" for kind in ("native", "adaptive")]
    pairs = [metrics.load(p / "fields.jsonl") for p in folders]
    a, b = pairs
    if [(r["stage"], r["progress"]) for r in a[0]] != [(r["stage"], r["progress"]) for r in b[0]]:
        raise ValueError("accepted stage or progress differs")
    model = observe.reference.Model(observe.read(root/"mesh.json"), observe.read(root/"ground.json"),
                                    observe.read(folders[0]/"case.json"))
    counts = [observe.partitions(p, len(model.seeds)) for p in folders]
    detail = observe.mismatch_rows(a, b, model.system, *counts)
    matrix = Counter(zip(a[1]["contact_status"].ravel(), b[1]["contact_status"].ravel(), strict=True))
    if sum(n for (i,j),n in matrix.items() if i != j) != len(detail):
        raise ValueError("strict categories and boundary diagnostics differ")
    ends = [i for i, r in enumerate(a[0]) if r["progress"] == 1.]
    result = {}
    for key in ("plastic_dissipation", "endpoint_work", "elastic_energy", "numerical_dissipation",
                "plastic_slip", "u", "traction", "water_load"):
        x, y = [fields[key] for _, fields in pairs]
        result[key] = dict(history=metrics.metric(x,y), endpoints=metrics.metric(x[ends],y[ends]))
        if key in ("plastic_dissipation", "endpoint_work"):
            dx, dy = [np.diff(v, axis=0, prepend=np.zeros_like(v[:1])) for v in (x,y)]
            result[key]["interval"] = metrics.metric(dx,dy)
    return dict(mechanical_steps=mechanical, water_steps=water, accepted=len(a[0]),
                confusion={f"{i}/{j}":n for (i,j),n in sorted(matrix.items())},
                strict_mismatches=len(detail), boundary=detail, quantities=result)


def measure(root):
    rows = [compare(root,m,w) for m,w in ((32,8),(32,16),(32,32),(16,32),(64,32))]
    with archive(SOURCE, "archive-index.json") as z:
        old = [json.loads(z.read(f"data/observation-final/mismatch-w{n}.json")) for n in (8,16,32)]
    return dict(mesh_tet4=12, contact_gp=6, categories=["open","stick","slip"],
                preserved_mismatch_counts=[len(v) for v in old], measured=rows,
                boundary_reclassification=False, physical_seconds=None, passed=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("output", type=Path)
    args=p.parse_args()
    value=measure(args.source)
    with args.output.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write("\n")
