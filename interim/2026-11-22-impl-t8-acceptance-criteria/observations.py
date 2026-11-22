"""Seal the small numerical comparison outputs and reproduce decisions from the seal."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile

from assess import assess
from sealed import HERE

FILES = (
    "a02/results.json", "a03/results.json", "a04/budget-summary.json",
    "t8-water/placement.json", "t8-water/affine.json",
    "b01/attached/measurements.json", "b01/mms/measurements.json",
    "b01/standard/measurements.json", "b01/modal/modes.json",
    "b02-patch/patch-results.json", "b02-area/comparison.json",
    "b03/controls.json", "b03/series.json", "b03/convergence.json",
    "c01-results.json", "c02-extended.json", "c04-results.json",
    "c03/results.json", "c03/local-replay/summary.json", "c03/matrix/restarts/verification.json",
    "c03/matrix/matrix-terminal.json", "c03/matrix/specifications.json",
    "c03/scope/scope.json", "c03/scope/mcc-pc105-open1/scope-extra.json",
    *(f"c04/{side}-m{m}-w{w}.time" for m,w in ((32,8),(32,16),(32,32),(16,32),(64,32))
      for side in ("native","adaptive")),
    "c01/native-128.time", "c01/native-256.time", "c01/native-512.time", "c03/matrix.time",
)


def pack(work, output):
    index = {}
    with zipfile.ZipFile(output / "observations.zip", "x") as z:
        for name in FILES:
            data = (work / name).read_bytes()
            json.loads(data)
            entry = zipfile.ZipInfo(name)
            z.writestr(entry, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
            index[name] = hashlib.sha256(data).hexdigest()
    for name in ("observations.zip", "acceptance.json", "comparisons.json"):
        index[name] = hashlib.sha256((output / name).read_bytes()).hexdigest()
    with (output / "observations.sha256").open("x", encoding="ascii", newline="\n") as stream:
        for name, sha in index.items():
            stream.write(f"{sha}  {name}\n")


def verify(source):
    lines = (source / "observations.sha256").read_text(encoding="ascii").splitlines()
    expected = {name: digest for digest,name in (line.split("  ",1) for line in lines)}
    if len(lines) != len(FILES)+3 or set(expected) != {*FILES, "observations.zip", "acceptance.json", "comparisons.json"}:
        raise ValueError("observation digest inventory differs")
    archive = source / "observations.zip"
    for name in ("observations.zip", "acceptance.json", "comparisons.json"):
        if hashlib.sha256((source/name).read_bytes()).hexdigest() != expected[name]:
            raise ValueError(f"changed comparison output {name}")
    with tempfile.TemporaryDirectory() as tmp, zipfile.ZipFile(archive) as z:
        if sorted(z.namelist()) != sorted(FILES):
            raise ValueError("unexpected comparison member")
        for name in FILES:
            data = z.read(name)
            if hashlib.sha256(data).hexdigest() != expected[name]:
                raise ValueError(f"changed comparison member {name}")
            path = Path(tmp) / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        decision, _ = assess(Path(tmp))
        if decision != json.loads((source / "acceptance.json").read_bytes()):
            raise ValueError("scope or outcome differs from measured fields")
    return len(FILES)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,default=HERE)
    p.add_argument("--pack",type=Path)
    args=p.parse_args()
    if args.pack:
        pack(args.pack,args.source)
    print(verify(args.source))
