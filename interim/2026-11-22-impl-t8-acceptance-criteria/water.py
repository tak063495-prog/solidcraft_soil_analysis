"""Independent cut-volume and saturated affine-head controls on the small columns."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import subprocess

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location(
    "placements", HERE.parent / "2026-11-20-impl-t2-soil-birth-hydraulic-init/study.py"
)
placement = importlib.util.module_from_spec(spec)
spec.loader.exec_module(placement)
controls = placement.t1


def solve(folder, solver, command, mesh, ground, case):
    folder.mkdir(parents=True, exist_ok=False)
    paths = [folder / f"{name}.json" for name in ("mesh", "ground", "case", "run")]
    for path, value in zip(paths[:3], (mesh, ground, case), strict=True):
        controls.write(path, value)
    with (folder / "solve.log").open("x") as log:
        subprocess.run([solver, command, *paths], stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=86400)
    run = controls.read(paths[-1])
    if not run["complete"]:
        raise ValueError(run.get("failure_reason"))
    return run


def placement_sweep(out, solver):
    rows = []
    for level in (0.625, 0.75, 0.875):
        mesh, ground, case = placement.prescribed("two_region")
        case["stages"][0]["hydraulic_initialization"]["level_m"] = level
        ground["initialization"]["water_level_m"] = level
        run = solve(out / f"placement-{level}", solver, "solve-stages", mesh, ground, case)
        event = next(e for e in run["transition_audits"] if e.get("placed_water"))
        n = 0.7 / 1.7
        lower, upper = level - 0.5, 1.0 - level
        expected = dict(lower_volume_m3=lower, upper_volume_m3=upper,
                        lower_water_m3=lower*n, upper_water_m3=upper*n*0.25,
                        lower_weight_kn=lower*((1-n)*27+n*10),
                        upper_weight_kn=upper*((1-n)*27+n*0.25*10))
        actual = {key: math.fsum(c[key] for c in event["placed_water"]) for key in expected}
        rows.append(dict(level_m=level, expected=expected, measured=actual,
                         difference={k: actual[k]-v for k, v in expected.items()},
                         zero_release_relative_force_residual=event["zero_release_relative_force_residual"],
                         whole_domain_budget=None))
    return rows


def affine_column(out, solver, slope, dt):
    case = dict(controls.cases())["continuous"]
    flow = case["flow"]
    for material in flow["materials"]:
        material["volume_source_per_s"] = 0.
    head = dict(head_m=2., reference_point_m=[0., 0., 0.], gradient=[0., 0., slope])
    flow["initial_head"] = head
    flow["analysis"]["periods"][0]["boundaries"] = [
        dict(face_set=side, condition=dict(kind="head", value=head)) for side in ("bottom", "top")
    ]
    flow["settings"].update(initial_time_step_s=dt, maximum_time_step_s=dt)
    mesh = controls.read(controls.EX / "joint-water-column-mesh.json")
    ground = controls.read(controls.EX / "ground-joint-water-mc.json")
    ground["initialization"].update(water_level_m=None)
    ground["initialization"]["geostatic"]["hydrostatic_face_sets"] = []
    folder = out / f"affine-{slope}-{dt}"
    run = solve(folder, solver, "solve-one-way", mesh, ground, case)
    budget = {}
    controls.record("continuous", mesh, run, folder, budget)
    records = []
    coords = np.array(mesh["nodes"])
    for frame in run["frames"]:
        if frame["phase"] != "transient":
            continue
        state = frame["state"]
        seconds = state["hydraulic"]["time_s"]
        traces = state["fluid"]["joint_traces"]
        forces = state["mechanical"]["joint_water"]["face_nodal_force_a_kn"]
        for t in traces:
            a, b = t["side_a"], t["side_b"]
            if seconds != t["time_s"] or a["face"] != b["face"] or a["owner"] == b["owner"]:
                raise ValueError("unpaired time or continuous face")
            z = t["elevation_m"]
            exact_head = 2. + slope*z
            pressure = 10.*(exact_head-z)
            face = t["face_index"]
            xyz = coords[mesh["face_sets"]["Seam:a"][face]]
            area = np.linalg.norm(np.cross(xyz[1]-xyz[0], xyz[2]-xyz[0]))/2
            cell_z = coords[mesh["tetrahedra"][a["owner"]], 2].mean()
            outward = math.copysign(1., z-cell_z)
            q = -1e-6*slope*area*outward
            expected_force = np.array([0., 0., -outward*pressure*area/3])
            records.append(dict(time_s=seconds, face=face, head_reference_m=exact_head,
                                side_a_head_error_m=a["total_head_m"]-exact_head,
                                side_b_head_error_m=b["total_head_m"]-exact_head,
                                pressure_error_kpa=t["applied_pressure_kpa"]-pressure,
                                expected_a_outflow_m3_s=q, a_flow_error_m3_s=a["outward_flow_m3_s"]-q,
                                b_flow_error_m3_s=b["outward_flow_m3_s"]+q,
                                max_nodal_force_error_kn=float(np.max(np.abs(np.array(forces[face])-expected_force)))))
    if {r["time_s"] for r in records} != {i*dt for i in range(1,round(.25/dt)+1)}:
        raise ValueError("incomplete physical-time series")
    return dict(slope=slope, physical_dt_s=dt, traces=records, water=budget["continuous"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("solver", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    solver = args.solver.resolve()
    controls.write(args.output / "placement.json", placement_sweep(args.output, solver))
    results = [affine_column(args.output, solver, slope, dt)
               for slope in (0., -0.2) for dt in (.125, .0625, .03125)]
    controls.write(args.output / "affine.json", results)
