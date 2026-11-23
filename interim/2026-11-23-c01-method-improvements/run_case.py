"""Execute one explicit method combination with the original event4 path budget."""

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
from methods import Methods, Seed, Sphere, STAGES


def context(steps):
    p = HERE.parent / "2026-12-12-c01-switch-event-localization/study.py"
    spec = importlib.util.spec_from_file_location("c01_method_study", p)
    study = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(study)
    mesh = study.read(ROOT / "examples/joint-water-column-tet10-mesh.json")
    ground = study.fem.normalize_ground(study.read(ROOT / "examples/ground-joint-water-mc.json"))
    stages = study.fem.normalize_stage_case(study.read(ROOT / "examples/stages-joint-water-cycle.json"))
    case = study.prior.configuration(stages, steps, 5e-5)
    case["settings"]["minimum_increment"] = min(
        case["settings"]["minimum_increment"], 1 / (steps * 65536)
    )
    control = dict(study.PATH_ERROR, event_proposal=True, minimum_contact_substeps=4)
    case["joint_system"]["global_step_control"] = dict(path_error=control)
    return study, mesh, ground, case, control


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def seal(folder, side, n):
    result = json.loads((folder / "result.json").read_bytes())
    if not result["complete"] or any(r["last_accepted"] != 1 for r in result["stages"]):
        raise ValueError("a complete coarse history is required")
    files = {
        "mesh": ROOT / "examples/joint-water-column-tet10-mesh.json",
        "ground": ROOT / "examples/ground-joint-water-mc.json",
        "case": folder / "stages.json",
        "result": folder / "result.json",
        "seed": folder / "seed.jsonl",
    }
    return dict(solver=side, steps=n, sha256={key: sha(path) for key, path in files.items()})


def execute(args):
    study, mesh, ground, case, control = context(args.steps)
    folder = args.folder
    folder.mkdir(parents=True)
    save = folder / "seed.jsonl"
    case_file = folder / "stages.json"
    seed = None
    if args.seed:
        source = args.seed
        expected = seal(source, args.solver, args.steps // 2)
        if (
            json.loads((source / "identity.json").read_bytes()) != expected
            or json.loads((source / "stages.json").read_bytes()) != context(args.steps // 2)[3]
        ):
            raise ValueError("seed solver, input or coarse time grid differs")
        seed = source / "seed.jsonl"
    if args.active:
        control["semismooth_contact"] = True
    case_file.write_text(json.dumps(case, allow_nan=False), encoding="utf-8")
    output = folder / "result.json"
    if args.solver == "native":
        options = []
        if seed or args.arc or args.active:
            opt = folder / "methods.json"
            opt.write_text(
                json.dumps(
                    dict(
                        seed=str(seed.resolve()) if seed else None,
                        arc_scale_m=0.02 if args.arc else None,
                    )
                ),
                encoding="utf-8",
            )
            options = ["--methods", opt]
        subprocess.run(
            [
                args.binary.resolve(),
                ROOT / "examples/joint-water-column-tet10-mesh.json",
                ROOT / "examples/ground-joint-water-mc.json",
                case_file,
                output,
                *([save] if args.save_seed else []),
                *options,
            ],
            check=True,
        )
        if args.save_seed:
            (folder / "identity.json").write_text(
                json.dumps(seal(folder, "native", args.steps)), encoding="utf-8"
            )
        return
    groups = {
        name: dict(
            stage=name, attempts=0, accepted=0, last_accepted=None, rejections=Counter(), max_force=0.0
        )
        for name in STAGES
    }

    def record(a):
        g = groups[a["stage"]]
        g["attempts"] += 1
        assert g["attempts"] <= 2000
        if a["accepted"]:
            g["accepted"] += 1
            g["last_accepted"] = a["progress_to"]
            audit = a["joint_path_error"]
            g["max_force"] = max(
                g["max_force"],
                a["relative_force_residual"],
                audit["midpoint_relative_force_residual"],
                audit["fine_relative_force_residual"],
            )
            assert (
                g["max_force"] <= 1e-10
                and audit["maximum_slip_error_m"] <= audit["slip_budget_m"]
                and audit["maximum_dissipation_error_kpa_m"] <= audit["dissipation_budget_kpa_m"]
            )
        else:
            g["rejections"][
                a.get(
                    "solver_failure",
                    "arc_correction"
                    if a.get("joint_path_error", {}).get("arc_correction")
                    else "path_error",
                )
            ] += 1

    sink = save.open("x", encoding="utf-8") if args.save_seed else None
    methods = Methods(
        sink,
        Seed(seed, len(mesh["nodes"])) if seed else None,
        args.active,
        Sphere(0.02) if args.arc else None,
    )
    try:
        study.reference._reference(
            mesh,
            ground,
            case,
            args.steps,
            1,
            7,
            False,
            True,
            maximum_m=5e-5,
            path_error=control,
            attempt_sink=record,
            methods=methods,
        )
        complete = True
        stop = None
    except AssertionError as e:
        if str(e) != "independent stage attempt budget":
            raise
        complete = False
        stop = str(e)
    finally:
        if sink:
            sink.close()
    output.write_text(
        json.dumps(
            dict(
                complete=complete,
                stop=stop,
                stages=list(groups.values()),
                seed_queries=methods.guesses,
                active_tangents=methods.tangents,
                arc_reviews=methods.arc.calls if methods.arc else 0,
            ),
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )

    if args.save_seed:
        (folder / "identity.json").write_text(
            json.dumps(seal(folder, "independent", args.steps)), encoding="utf-8"
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    p.add_argument("--steps", type=int, choices=(256, 512, 1024), required=True)
    p.add_argument("--solver", choices=("native", "independent"), required=True)
    p.add_argument("--seed", type=Path)
    p.add_argument("--arc", action="store_true")
    p.add_argument("--active", action="store_true")
    p.add_argument("--save-seed", action="store_true")
    p.add_argument("--binary", type=Path, default=Path("target/release/examples/c01_methods"))
    execute(p.parse_args())
