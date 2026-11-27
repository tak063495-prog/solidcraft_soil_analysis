"""Collect actual elapsed time, process RSS and limits from completed runs."""
import hashlib
import json
from pathlib import Path
import sys
from plans import HERE,ROUTES,save


def profile(path):
    d=json.loads(path.read_text().splitlines()[-1])
    assert d["exit_code"]==0 and d["seconds"]<86400
    return d


def run(out):
    files=[out/"checks"/f"face-{s}.log" for s in ("plus","minus")]
    files.append(out/"checks/local-lu.log")
    r=out/"c01/reload"
    files.extend([r/"prefix-resource.log",r/"resource.log"])
    for route in ROUTES:
        files.extend(out/"c01/alternate"/route/f"{n}-resource.log" for n in (256,512,1024))
    logs={str(f.relative_to(out)):profile(f) for f in files}
    prefix_and_reload=sum(logs[str(f.relative_to(out))]["seconds"] for f in (r/"prefix-resource.log",r/"resource.log"))
    assert prefix_and_reload<86400
    memory={k:int(Path("/sys/fs/cgroup/memory."+k).read_text()) for k in ("max","peak")}
    assert memory["peak"]<=memory["max"]<=40_000_000_000
    inputs=[HERE/"instrument.py",HERE/"face_chain.py",HERE/"local_lu.rs",HERE/"local_lu.py",
            *(out/"c01/alternate"/route/"plan.json" for route in ROUTES)]
    save(out/"checks/resources.json",dict(mesh_elements=dict(joint=12,nc=96),stage_cap=1787,
         proportional_caps={"12":2065,"96":16520},default_cap=2000,profile=logs,
         prefix_and_reload_seconds=prefix_and_reload,cgroup_memory_bytes=memory,
         sha256={str(f.relative_to(out) if f.is_relative_to(out) else f.name):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs}))


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve())
