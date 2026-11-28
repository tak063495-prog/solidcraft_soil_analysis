"""Build three 12 Tet10 routes from common inputs for the one-axis comparisons."""
import argparse
import json
from pathlib import Path

from c03 import Model, read, write
from c03diag_pc import export_native
from c03diag_reference import run
from c03diag_study import binary, configuration, initialize, invoke, schedule

SIZES = {"m0c1": (1e-4, 5e-5), "m1c0": (5e-5, 1e-4), "m1c1": (5e-5, 5e-5)}


def solve(root, label, solver, *, resolve_edges=False):
    folder = root / f"{label}-{solver}"
    folder.mkdir()
    material, contact = SIZES[label]
    mesh, ground = [read(root / f"{name}.json") for name in ("mesh", "ground")]
    if len(mesh.get("tetrahedra10", [])) != 12:
        raise ValueError("12 Tet10 required")
    case = configuration(root, 32, material, contact)
    write(folder/"stages.json", case)
    if solver == "native":
        process = invoke([binary(root), root/"mesh.json", root/"ground.json", folder/"stages.json",
                          folder/"raw"], folder/"solve.log", 14400)
        write(folder/"process.json", process)
        if process["exit_code"]:
            raise RuntimeError("native route incomplete")
        write(folder/"schedule.json", schedule(read(folder/"raw/run.json"), 32))
        export_native(root, folder, 32)
        return
    model = Model(mesh, ground, case)
    with (folder/"fields.jsonl").open("x", encoding="utf-8") as anchors, (
            folder/"substeps.jsonl").open("x", encoding="utf-8") as extras:
        def accepted(row, fields):
            stream = anchors if row["progress"] * 32 == round(row["progress"] * 32) else extras
            print(json.dumps(dict(row=row, fields={k:v.tolist() for k,v in fields.items()}),
                             allow_nan=False), file=stream, flush=True)
        result = run(model, 32, material, contact, accepted, cutback=True, resolve_edges=resolve_edges)
    write(folder/"result.json", result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--initialize", type=Path, metavar="NATIVE_PROBE")
    parser.add_argument("--route", choices=SIZES)
    parser.add_argument("--solver", choices=("native", "independent"))
    parser.add_argument("--resolve-edges", action="store_true",
                        help="converge coactive independent MC shear planes to their classification width")
    args = parser.parse_args()
    if args.resolve_edges and args.solver != "independent":
        parser.error("--resolve-edges requires an independent MC route")
    if args.initialize and not (args.route or args.solver):
        initialize(args.root, args.initialize.resolve(), "tet10", 30., open_factor=.12)
    elif args.route and args.solver and not args.initialize:
        solve(args.root, args.route, args.solver, resolve_edges=args.resolve_edges)
    else:
        parser.error("select --initialize or both --route and --solver")
