"""Run the P2 endpoint solver on a preselected asymmetric reversal cycle."""
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
P2 = HERE.parent / "2026-11-25-remaining5-p2"
ROOT = HERE.parents[1]
sys.path.insert(0, str(P2))
from c01_native import Fields, context


def run(out, n, fixed_budget):
    out.mkdir(parents=True)
    _, _, _, case, control = context(n)
    case["settings"]["attempt_budget_policy"] = "mesh_proportional"
    control.update(semismooth_contact=True, contact_history_bisection=True,
                   localize_friction_events=True, resolve_slip_path=True,
                   reequilibrate_endpoints=True)
    if fixed_budget:
        control["slip_scale_m"] = n**2 * (1e-4 / 1024**2)
    targets = {"reverse": -0.075, "reload": 0.06}
    for stage in case["stages"]:
        if stage["name"] in targets:
            stage["target_factor"] = targets[stage["name"]]
    (out / "stages.json").write_text(json.dumps(case, indent=2) + "\n")
    (out / "methods.json").write_text(json.dumps(dict(
        full_field_grid=n, dissipation_history=str((out / "dissipation.jsonl").resolve()))))
    args = [ROOT / "target/release/examples/c01_methods",
            ROOT / "examples/joint-water-column-tet10-mesh.json",
            ROOT / "examples/ground-joint-water-mc.json", out / "stages.json",
            out / "result.json", "--methods", out / "methods.json"]
    with Fields(out / "fields", n, "native") as fields, subprocess.Popen(args, stdout=subprocess.PIPE) as proc:
        try:
            for line in proc.stdout:
                fields.native(*json.loads(line))
        except BaseException:
            proc.kill()
            raise
        if proc.wait():
            raise subprocess.CalledProcessError(proc.returncode, args)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("n", type=int, choices=[256, 512, 1024])
    parser.add_argument("--fixed-budget", action="store_true")
    args = parser.parse_args()
    run(args.out, args.n, args.fixed_budget)
