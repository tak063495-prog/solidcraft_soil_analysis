"""Collect path costs and identify prior finite-element inputs by content."""
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def collect(d):
    raw={}
    for file in sorted(d.rglob("*-resource.log")):
        if file.stat().st_size:
            raw[str(file.relative_to(d))]=json.loads(file.read_text().splitlines()[-1])
    wall=lambda key:raw[key]["seconds"]
    prefix=wall("c01/reload/prefix-resource.log")
    seconds={f"reverse+{name}":prefix+wall(f"c01/reload/{name}-resource.log") for name in ("observation","narrow")}
    seconds.update({f"alternate{n}":wall(f"c01/alternate/{n}-resource.log") for n in (256,512,1024)})
    inputs=[ROOT/"examples"/(n+".json") for n in ("joint-water-column-tet10-mesh","ground-joint-water-mc","stages-joint-water-cycle")]
    for sign in ("plus","minus"):
        base=HERE.parent/"2026-11-26-remaining-p8/c02"/sign
        inputs.extend(base/f for f in ("kinematics.npz","slip.npz","tiles.json","residual.npz"))
    inputs.extend(HERE.parent/f"2026-11-26-remaining-p7/c01/{i}/input.json" for i in ("original1024","alternate256","alternate512","alternate1024"))
    inputs.extend((HERE.parent/"2026-11-26-remaining-p9/c01/reload/reverse-end.json.gz",HERE.parent/"2026-11-26-remaining-p8/nc/parent-632.json.gz"))
    sha={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs}
    cg=Path("/sys/fs/cgroup")
    limit=int((cg/"memory.max").read_text())
    peak=int((cg/"memory.peak").read_text())
    assert 0<peak<=limit<=40_000_000_000 and max(seconds.values())<=86400
    out=dict(base="6a6b5f83887f5f776ccc03b0e3c63bbb4caf59b0",cgroup_hard_limit_bytes=limit,cgroup_peak_bytes_at_audit=peak,
             mesh_elements=dict(joint=12,nc=96),joint_attempt_cap=1787,mesh_budget=dict(joint=2065,nc=16520),
             default_cap=2000,route_seconds=seconds,profile_logs=raw,inputs_sha256=sha)
    (d/"checks/resources.json").write_text(json.dumps(out,indent=2)+"\n")


if __name__=="__main__":
    collect(Path(sys.argv[1]).resolve())
