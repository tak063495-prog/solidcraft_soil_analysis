"""Collect route measurements and fingerprints of the frozen sources."""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
P13 = HERE.with_name("2026-11-27-remaining-p13")
P14 = HERE.with_name("2026-11-27-remaining-p14")
P15 = HERE.with_name("2026-11-28-remaining-p15")


def record(path):
    item = json.loads(path.read_text().splitlines()[-1])
    assert item["exit_code"]==0 and item["seconds"]<86400 and item["peak_rss_kib"]*1024<40_000_000_000
    return item


def digest(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def run(work,out):
    checks=out/"checks"
    runs = {k:record(checks/(k+".log")) for k in
          ["nc","face-plus","face-minus",*[f"{s}-{n}" for s in ("outer","middle") for n in (256,512,1024)]]}
    for k,s in (("reload-prefix","prefix-resource.log"),("reload","resource.log")):
        runs[k]=record(out/"c01/reload"/s)
    seconds = runs["reload-prefix"]["seconds"]+runs["reload"]["seconds"]
    assert seconds<86400
    memory_cap = int(Path("/sys/fs/cgroup/memory.max").read_text())
    assert memory_cap<=40_000_000_000
    files = [P13/"nc/global.jsonl",P15/"nc/schur.json",P15/"c01/reload/trace.jsonl.gz",
             P15/"c01/alternate/common-middle.jsonl.gz",P15/"c02/plus.json",P15/"c02/minus.json"]
    files += [P13/f"nc/{s}-end.json.gz" for s in ("ordinary","predicted")]
    for name,p in (("outer",P15),("middle",P14)):
        for n in (256,512,1024):
            files += [p/f"c01/alternate/{name}/{n}/{s}.json" for s in ("input","manifest")]
    bins = work/"probe/target/release/examples"
    data = dict(start_sha="bdbfaa22900e14310e6589251a6445caccdc29dc",approval_jst="2026-11-28 08:39",
        runs=runs,reload_including_prefix_seconds=seconds,wall_cap_seconds=86400,
        cgroup_limit_bytes=memory_cap,cgroup_peak_bytes=int(Path("/sys/fs/cgroup/memory.peak").read_text()),
        contact_elements=12,nc_elements=96,nc_active_elements=48,controlled_stage_limit=1787,
        mesh_budget_12=2065,mesh_budget_96=16520,fixed_default=2000,
        input_sha256={str(p.relative_to(HERE.parents[1])):digest(p) for p in files},
        probe_sha256={s:digest(bins/s) for s in ("c01_methods","local_lu","p5_homogeneous","p5_point")})
    (checks/"resources.json").write_text(json.dumps(data,indent=2)+"\n")


if __name__=="__main__":
    run(*map(Path,sys.argv[1:]))
