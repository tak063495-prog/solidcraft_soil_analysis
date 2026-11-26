"""Measure charged path runtimes, input identity and the shared hard memory ceiling."""
import hashlib
import json
from pathlib import Path
import sys


def main(out):
    names=[f"c01/alternate/{n}-resource.log" for n in (256,512,1024)]
    names+=["c01/reload/prefix-resource.log","c01/reload/entrance-resource.log","checks/nc.log",
            "checks/face-plus.log","checks/face-minus.log"]
    profiles={f:json.loads((out/f).read_text().splitlines()[-1]) for f in names}
    assert all(v["exit_code"]==0 for v in profiles.values())
    seconds={f:v["seconds"] for f,v in profiles.items() if "reload" not in f}
    seconds["reload_prefix_and_entrance"]=sum(profiles[f]["seconds"] for f in names if "reload" in f)
    stops=[json.loads(f.read_text().splitlines()[-1]) for f in out.glob("checks/nc-*.log")]
    seconds["nc_all_invocations"]=seconds["checks/nc.log"]+sum(r["seconds"] for r in stops)
    limit=int(Path("/sys/fs/cgroup/memory.max").read_text())
    peak=int(Path("/sys/fs/cgroup/memory.peak").read_text())
    assert max(seconds.values())<86400 and peak<=limit<=40_000_000_000
    files=["c01/reload/plan.json","c01/alternate/plan.json","nc/global.jsonl"]
    summary=dict(mesh_elements=dict(joint=12,nc=96),mesh_caps=dict(joint=2065,nc=16520),
                 joint_attempt_cap=1787,default=2000,cgroup_limit_bytes=limit,cgroup_peak_bytes=peak,
                 runtime_seconds=seconds,profiles=profiles,
                 sha256={f:hashlib.sha256((out/f).read_bytes()).hexdigest() for f in files})
    (out/"checks/resources.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__=="__main__":
    main(Path(sys.argv[1]))
