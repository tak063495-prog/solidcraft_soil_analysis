"""Summarize per-path profiling, charged trials and local verification logs."""
import hashlib
import json
from pathlib import Path
import re
import sys

HERE=Path(__file__).resolve().parent
BASE="dd57540d627fea2b1693023be591d73451fe4b59"


def run(work,out):
    checks=out/"checks"
    logs={s:checks/(s+".log") for s in ("face-plus","face-minus","nc",*[
        f"{r}-{n}" for r in ("middle","outer") for n in (256,512,1024)])}
    for key,file in (("reload","resource.log"),("reload-prefix","prefix-resource.log")):
        logs[key]=out/"c01/reload"/file
    measured={k:json.loads(p.read_text().splitlines()[-1]) for k,p in logs.items()}
    assert all(v["exit_code"]==0 and v["seconds"]<86400 for v in measured.values())
    reload_seconds=sum(measured[k]["seconds"] for k in ("reload","reload-prefix"))
    assert reload_seconds<86400
    root=Path("/sys/fs/cgroup")
    def digest(path):
        with path.open("rb") as f:
            return hashlib.file_digest(f,"sha256").hexdigest()
    usage=dict(base_sha=BASE,approval_jst="2026-11-28 13:27",cgroup_limit_bytes=int((root/"memory.max").read_text()),
               cgroup_peak_bytes=int((root/"memory.peak").read_text()),wall_limit_seconds=86400,
               runs=measured,reload_with_prefix_seconds=reload_seconds,contact_elements=12,nc_tet10=96,
               nc_active=48,charged_stage_limit=1787,budget_at_12=2065,default_limit=2000,
               probe_sha256={k:digest(work/"probe/target/release/examples"/k)
                   for k in ("c01_methods","local_lu","p5_homogeneous","p5_point")})
    assert usage["cgroup_limit_bytes"]<=40_000_000_000
    (checks/"resources.json").write_text(json.dumps(usage,indent=2)+"\n")
    rust_log=(checks/"rust.log").read_text()
    counts=dict(rust=sum(int(n) for n in re.findall(r"test result: ok\. (\d+) passed",rust_log)))
    assert counts["rust"]>0 and "FAILED" not in rust_log
    for n in range(9,18):
        log=(checks/f"p{n}-tests.log").read_text()
        assert log.rstrip().endswith("OK")
        counts[f"p{n}"]=int(re.search(r"Ran (\d+) tests",log)[1])
    (checks/"test-counts.json").write_text(json.dumps(dict(suites=counts,total=sum(counts.values())),indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
