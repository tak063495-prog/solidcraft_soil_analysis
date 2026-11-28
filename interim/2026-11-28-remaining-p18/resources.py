"""Summarize exit codes, charged runs, test counts and the observed cgroup ceiling."""
import hashlib
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent


def read_json(path):
    return json.loads(path.read_text())


def save(path, value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


def record(work,out):
    checks = out/"checks"
    profiles = {s:checks/(s+".log") for s in ("face-minus","face-plus","face-minus-328","nc",
                    "outer-256","outer-512","outer-1024","middle-256","middle-512","middle-1024")}
    profiles.update(prefix=out/"c01/reload/prefix-resource.log",shadow=out/"c01/reload/resource.log",
                    control=out/"c01/reload/control/resource.log")
    metrics = {s:json.loads(p.read_text().splitlines()[-1]) for s,p in profiles.items()}
    assert all(v["exit_code"]==0 and v["seconds"]<86400 for v in metrics.values())
    whole = {k:metrics["prefix"]["seconds"]+metrics[k]["seconds"] for k in ("shadow","control")}
    assert max(whole.values())<86400
    cgroup = Path("/sys/fs/cgroup")
    cap = int((cgroup/"memory.max").read_text())
    peak = int((cgroup/"memory.peak").read_text())
    events = {k:int(v) for k,v in (line.split() for line in (cgroup/"memory.events").read_text().splitlines())}
    assert peak<=cap<=40_000_000_000
    save(checks/"resources.json",dict(base_sha=(checks/"base.txt").read_text().strip(),
        approval_jst="2026-11-28 16:51",memory_limit_bytes=cap,cgroup_peak_bytes=peak,
        cgroup_memory_events=events,time_limit_seconds=86400,mesh_elements=dict(contact=12,nc=96),runs=metrics,
        original_including_prefix_seconds=whole,probe_sha256={
            k:hashlib.sha256((work/"probe/target/release/examples"/k).read_bytes()).hexdigest()
            for k in ("c01_methods","local_lu","p5_homogeneous","p5_point")}))
    suites = {}
    for i in range(9,19):
        log=(checks/f"p{i}-tests.log").read_text()
        m=re.search(r"Ran (\d+) tests? in ",log)
        assert m and log.rstrip().endswith("OK")
        suites[f"p{i}"]=int(m[1])
    log=(checks/"rust.log").read_text()
    rust=re.findall(r"test result: ok\. (\d+) passed; 0 failed; 0 ignored; 0 measured; 0 filtered out",log)
    assert len(rust)==3
    save(checks/"test-counts.json",dict(rust=sum(map(int,rust)),python=suites,
        total=sum(map(int,rust))+sum(suites.values())))


if __name__ == "__main__":
    record(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
