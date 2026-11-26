"""Record route budgets, locked inputs and the active cgroup ceiling."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
P8=HERE.parent/"2026-11-26-remaining-p8"
P7=HERE.parent/"2026-11-26-remaining-p7"


def main(out):
    logs=list((out/"checks").glob("*.log"))+list((out/"c01").rglob("*-resource.log"))
    jobs={}
    for path in sorted(logs):
        records=[json.loads(l) for l in path.read_text().splitlines() if l.startswith('{"exit_code":')]
        if records: jobs[str(path.relative_to(out))]=records
    inputs=[P8/"c01/reverse-end-checkpoint.json.gz",P8/"c01/predictor.json",P8/"nc/parent-457.json.gz",
            P8/"nc/parent-632.json.gz",P7/"nc/continuation.jsonl",P8/"c02/comparison.json"]
    for sign in ("plus","minus"):
        inputs.extend((P8/"c02"/sign).glob("*.npz"))
        inputs.append(P8/"c02"/sign/"tiles.json")
    for n in (256,512,1024):
        inputs.extend([P7/f"c01/alternate{n}/input.json",P7/f"c01/alternate{n}/trace.jsonl.gz"])
    inputs.extend([P7/"c01/original1024/input.json",ROOT/"examples/ground-joint-water-mc.json",
                   ROOT/"examples/joint-water-column-tet10-mesh.json",ROOT/"examples/stages-joint-water-cycle.json"])
    prefix=jobs["checks/prefix-resource.log"][-1]["seconds"]
    paths={name:prefix+jobs[f"c01/reload/{name}-resource.log"][-1]["seconds"] for name in ("baseline","predicted","windows")}
    paths.update({f"alternate{n}":jobs[f"c01/alternate/{n}-resource.log"][-1]["seconds"] for n in (256,512,1024)})
    paths.update({f"matched{n}":jobs[f"c01/alternate/matched-{n}-resource.log"][-1]["seconds"] for n in (256,512,1024)})
    info=dict(start_sha="0e26bda4c30cb1125a2bdf386bc655c86376a370",
              preregistration_sha="af32241311ba20403522e10a8ed8f8f8ad6894ab",
              approval_jst="2026-11-26 18:04",python=sys.version.split()[0],numpy=np.__version__,
              rust=subprocess.check_output(["rustc","--version"],text=True).strip(),
              cgroup_limit_bytes=int(Path("/sys/fs/cgroup/memory.max").read_text()),
              cgroup_peak_bytes=int(Path("/sys/fs/cgroup/memory.peak").read_text()),
              joint_path_seconds=paths,jobs=jobs,
              input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    assert max(paths.values())<86400 and info["cgroup_limit_bytes"]<=40_000_000_000
    (out/"checks/resources.json").write_text(json.dumps(info,indent=2)+"\n")


if __name__=="__main__":
    main(Path(sys.argv[1]).resolve())
