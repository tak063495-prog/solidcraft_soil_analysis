"""Replay the common prefix of both archived fixed-budget cases from geostatic."""
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
P3 = HERE.parent / "2026-11-25-remaining-p3/c01/fixed"
P4 = HERE.parent / "2026-11-25-remaining-p4/c01"


def run(out, binary, n):
    out.mkdir(parents=True)
    original = P4 / f"original{n}-input.json"
    alternate = P3 / f"case-{n}.json"
    cases = [json.loads(p.read_bytes()) for p in (original, alternate)]
    for case in cases:
        assert [s["name"] for s in case["stages"][:3]] == ["shear", "open", "recontact"]
        case["stages"] = case["stages"][:3]
    assert cases[0] == cases[1]
    input_file = out / "prefix.json"
    input_file.write_text(json.dumps(cases[0]))
    cmd = [binary, ROOT / "examples/joint-water-column-tet10-mesh.json",
           ROOT / "examples/ground-joint-water-mc.json", input_file, out / "result.json"]
    with gzip.open(out / "trials.jsonl.gz", "wb", compresslevel=4) as sink:
        with subprocess.Popen(cmd, stderr=subprocess.PIPE) as proc:
            for line in proc.stderr:
                sink.write(line)
            if proc.wait():
                raise RuntimeError(f"replay failed, see {sink.name}")
    result = json.loads((out / "result.json").read_bytes())
    previous = json.loads((P4 / f"original{n}-result.json").read_bytes())
    alt_previous = json.loads((P3 / "summary.json").read_bytes())["routes"][str(n)]
    assert result["complete"] and result["stages"][:4] == previous["stages"][:4]
    assert result["stages"][:4] == alt_previous["stages"][:4]
    with gzip.open(out / "trials.jsonl.gz", "rt") as stream:
        rows = list(map(json.loads, stream))
    extra = sum(row["extra_trial"] for row in rows)
    assert result["stages"][3]["attempts"] + extra <= 2065
    kept = [r for r in rows if r["accepted"] and (r["end"] * 256).is_integer()]
    np.testing.assert_array_equal([r["end"] * 256 for r in kept], np.arange(129, 257))
    slip = [r["after_f"]["slip"] for r in kept]
    for path in (P4 / "original-fixed-slip.npz", P3 / "slip-common.npz"):
        with np.load(path) as data:
            np.testing.assert_array_equal(slip, data[f"n{n}"][3*256+128:4*256])
    (out / "audit.json").write_text(json.dumps(dict(common_prefix=True,
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (original, alternate)}, attempts=result["stages"][3]["attempts"],
        additional_material_and_global_trials=extra, matched_anchors_per_route=len(kept),
        archived_prefix_stage_summaries_identical=True), indent=2)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]), Path(sys.argv[2]).resolve(), int(sys.argv[3]))
