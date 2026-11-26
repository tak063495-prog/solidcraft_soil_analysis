"""Summarize recorded wall times, observed charges and the runtime memory limit."""
import hashlib
import json
from pathlib import Path
import sys

from plans import save


def run(out):
    profiles={}
    for p in sorted(out.rglob("*.log")):
        for s in p.read_text(errors="replace").splitlines():
            if s.startswith('{"exit_code"'):
                profiles[str(p.relative_to(out))]=json.loads(s)
    reload_time=sum(profiles["c01/reload/"+f]["seconds"] for f in ("prefix-resource.log","resource.log"))
    path_times={"original_prefix_and_reload":reload_time}
    for key,v in profiles.items():
        if key.startswith("c01/alternate/"):
            path_times[key]=v["seconds"]
    path_times["both_nc_paths"]=profiles["checks/nc.log"]["seconds"]
    assert all(t<86400 for t in path_times.values())
    cgroup=Path("/sys/fs/cgroup")
    cap=int((cgroup/"memory.max").read_text())
    peak=int((cgroup/"memory.peak").read_text())
    assert peak<=cap<=40_000_000_000
    hashes={}
    for name in ("design.json","c01/reload/plan.json","c01/alternate/point14/plan.json",
                 "c01/alternate/point26/plan.json","nc/global.jsonl"):
        hashes[name]=hashlib.sha256((out/name).read_bytes()).hexdigest()
    save(out/"checks/resources.json",dict(mesh=dict(joint=12,nc_tet10=96),budget=dict(joint_stage=1787,
        B12=2065,B96=16520,default=2000),cgroup_limit_bytes=cap,cgroup_peak_bytes=peak,
        path_times_seconds=path_times,profiles=profiles,sha256=hashes))


if __name__=="__main__":
    run(Path(sys.argv[1]))
