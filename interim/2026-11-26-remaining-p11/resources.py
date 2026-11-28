"""Inventory completed routes at their actual mesh and capture the hard resource limit."""
import hashlib
import json
from pathlib import Path
import sys

def main(d):
    profiles={}
    for f in [*d.glob("c01/*/*-resource.log"),*d.glob("checks/face-*.log"),d/"checks/nc.log"]:
        rows=f.read_text().splitlines()
        r=json.loads(rows[-1])
        assert r["exit_code"]==0 and r["seconds"]<86400
        profiles[str(f.relative_to(d))]=r
    def seconds(name):
        return profiles[name]["seconds"]
    # Each resumed reload route includes the cost of its shared reverse prefix.
    times={s:seconds("c01/reload/prefix-resource.log")+seconds("c01/reload/"+s+"-resource.log")
           for s in ("cross","ordinary","half")}
    times.update({f"alternate-{n}":seconds(f"c01/alternate/{n}-resource.log") for n in (256,512,1024)})
    times.update({f"observe-{n}":seconds(f"c01/alternate/observe-{n}-resource.log") for n in (256,512,1024)})
    times["two_nc_paths_with_replay"]=seconds("checks/nc.log")
    limit=int(Path("/sys/fs/cgroup/memory.max").read_text())
    peak=int(Path("/sys/fs/cgroup/memory.peak").read_text())
    assert max(times.values())<86400 and peak<=limit<=40_000_000_000
    record=dict(root="5bb94fc6deb18696c06d08e188b53c766417ac66",
        mesh_elements=dict(joint=12,nc=96),mesh_budget=dict(joint=2065,nc=16520),joint_cap=1787,default_cap=2000,
        cgroup_limit_bytes=limit,cgroup_peak_bytes=peak,profiles=profiles,seconds=times,
        sha256={str(f.relative_to(d)):hashlib.sha256(f.read_bytes()).hexdigest() for f in
                (d/"c01/reload/plan.json",d/"c01/alternate/plan.json",d/"nc/global.jsonl")})
    (d/"checks/resources.json").write_text(json.dumps(record,indent=2)+"\n")

if __name__=="__main__":
    main(Path(sys.argv[1]))
