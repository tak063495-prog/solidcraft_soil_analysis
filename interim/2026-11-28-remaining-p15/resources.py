"""Record input bindings, route wall times and the enforced memory limit."""
import hashlib
import json
from pathlib import Path
import re
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
P14=HERE.with_name("2026-11-27-remaining-p14")
P13=HERE.with_name("2026-11-27-remaining-p13")


def main(out):
    profiles={}
    for pattern in ("checks/face-plus.log","checks/face-minus.log","c01/reload/*resource.log","c01/alternate/outer/*-resource.log"):
        for p in sorted(out.glob(pattern)):
            match=re.search(r'\{[^\n]*"peak_rss_kib"[^\n]*\}',p.read_text())
            if not match:
                raise ValueError(f"missing route profile: {p}")
            r=json.loads(match.group())
            assert r["exit_code"]==0 and r["seconds"]<86400
            profiles[str(p.relative_to(out))]=r
    reload_seconds=sum(r["seconds"] for k,r in profiles.items() if k.startswith("c01/reload/"))
    assert reload_seconds<86400
    sha={}
    inputs=[*sorted(HERE.glob("*.py")),*sorted(HERE.glob("*.sh")),
            P13/"nc/direction-check.json",P14/"nc/local-lu.json",P14/"face_chain.py",P13/"inverse.py",
            P14/"c01/alternate/outer/plan.json",P14/"c01/alternate/middle/plan.json",
            HERE.with_name("2026-11-26-remaining-p10")/"c01/reload/reverse-end.json.gz",
            HERE.with_name("2026-11-27-remaining-p12")/"c01/reload/plan.json"]
    for p in inputs:
        sha[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    cg=Path("/sys/fs/cgroup")
    mem={k:int((cg/("memory."+k)).read_text()) for k in ("peak","max")}
    assert mem["max"]<=40_000_000_000
    (out/"checks/resources.json").write_text(json.dumps(dict(base="2bbab82088846137c6384a48b926eb24229b2137",
        approved="2026-11-28T05:34:00+09:00",profiles=profiles,reload_including_prefix_seconds=reload_seconds,memory_cgroup_bytes=mem,
        stage_cap=1787,default=2000,proportional_cap={"12":2065,"96":16520},
        mesh_elements={"contact":12,"nc":96},sha256=sha),indent=2)+"\n")


if __name__=="__main__":
    main(Path(sys.argv[1]))
