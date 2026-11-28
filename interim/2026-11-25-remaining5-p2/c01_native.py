"""Replay the native P1 case with only endpoint re-equilibration selected."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0,str(HERE.parent/"2026-11-24-c01-event4-retry"))
from c01_retry_fields import Fields, context


def run(folder, n, binary):
    folder.mkdir(parents=True)
    _, _, _, case, control = context(n)
    case["settings"]["attempt_budget_policy"] = "mesh_proportional"
    control.update(semismooth_contact=True,contact_history_bisection=True,
                   localize_friction_events=True,resolve_slip_path=True,reequilibrate_endpoints=True)
    stages = folder/"stages.json"
    stages.write_text(json.dumps(case,indent=2)+"\n")
    options = folder/"methods.json"
    options.write_text(json.dumps(dict(full_field_grid=n,dissipation_history=str((folder/"dissipation.jsonl").resolve()))))
    command = [binary,ROOT/"examples/joint-water-column-tet10-mesh.json",
               ROOT/"examples/ground-joint-water-mc.json",stages,folder/"result.json","--methods",options]
    with Fields(folder/"fields",n,"native") as fields, subprocess.Popen(command,stdout=subprocess.PIPE) as proc:
        try:
            for line in proc.stdout:
                fields.native(*json.loads(line))
        except BaseException:
            proc.kill()
            raise
        if proc.wait():
            raise subprocess.CalledProcessError(proc.returncode,command)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder",type=Path)
    parser.add_argument("n",type=int,choices=[256,512,1024])
    parser.add_argument("binary",type=Path)
    args=parser.parse_args()
    run(args.folder,args.n,args.binary.resolve())
