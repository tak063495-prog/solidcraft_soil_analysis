"""Run the C-03 one-axis controls on the existing wet MCC pc105 material."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-23-impl-root-causes-c"))
from c03_cloud import SIZES
from c03diag_study import binary, configuration, initialize, read, write
from c03diag_mcc_reference import Model, run
from independent_work import AcceptedWork


def prepare(root, probe):
    base = root / "basis"
    initialize(base, probe.resolve(), "tet10", 30., open_factor=.12)
    ground = read(base / "ground.json")
    for material in ground["materials"]:
        del material["mohr_coulomb"]
        material["mcc"] = {"critical_state_m": 1.2, "lambda": .2, "kappa": .04,
                           "poisson": .3, "reference_pressure_kpa": 1.,
                           "specific_volume_ncl_at_reference": 3.5}
    for seed in ground["initialization"]["geostatic"]["seed"]["states"].values():
        for key in ("equivalent_plastic_strain", "plastic_work_kpa", "last_return"):
            del seed[key]
        seed["preconsolidation_pressure_kpa"] = 105.
    write(root / "ground.json", ground)
    for name, (material, contact) in SIZES.items():
        folder = root / name
        folder.mkdir()
        write(folder / "stages.json", configuration(base, 32, material, contact))


def solve(root, name, solver):
    import subprocess
    base = root / "basis"
    folder = root / name
    if solver == "native":
        subprocess.run([binary(base), base / "mesh.json", root / "ground.json", folder / "stages.json",
                        folder / "raw", "--diagnose"], check=True, timeout=14400)
        return
    model = Model(read(base / "mesh.json"), read(root / "ground.json"), read(folder / "stages.json"))
    work = AcceptedWork(model)
    with (folder / "independent-work.jsonl").open("x", encoding="utf-8") as stream:
        def accepted(row, fields):
            entry = work.observe(row, fields)
            entry["relative_force_residual"] = row["relative_force_residual"]
            print(json.dumps(entry, allow_nan=False), file=stream, flush=True)
        result = run(model, 32, *SIZES[name], accepted)
    if not result["complete"] or result["reference_uses_native_results"]:
        raise ValueError("independent MCC comparison incomplete")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("prepare").add_argument("probe", type=Path)
    run_case = sub.add_parser("solve")
    run_case.add_argument("name", choices=SIZES)
    run_case.add_argument("solver", choices=("native", "independent"))
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.root, args.probe)
    else:
        solve(args.root, args.name, args.solver)
