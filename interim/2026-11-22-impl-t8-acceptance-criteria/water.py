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
AFFINE_SLOPES = (0., -0.2)
AFFINE_STEPS = (.125, .0625, .03125)
TRI3_POINTS = range(3)


def affine_fixture():
    mesh = controls.read(controls.EX / "joint-water-column-mesh.json")
    return mesh, dict(controls.cases())["continuous"]


def affine_faces(mesh, case):
    joint, = case["mechanics"]["joint_system"]["joints"]
    faces = mesh["face_sets"][joint["side_a"]]
    if (not faces or mesh.get("tetrahedra10") or any(len(face) != 3 for face in faces)
            or len(faces) != len(mesh["face_sets"][joint["side_b"]])):
        raise ValueError("affine comparison requires paired Tri3 faces")
    # Tet4 JOINTs use the three-point Tri3 rule in elements::face_points_order.
    return joint["id"], faces


def affine_times(case, dt):
    period, = case["flow"]["analysis"]["periods"]
    return {i * dt for i in range(1, round(period["end_time_s"] / dt) + 1)}


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
    mesh, case = affine_fixture()
    joint_id, faces = affine_faces(mesh, case)
    expected_points = {(joint_id, face, point) for face in range(len(faces)) for point in TRI3_POINTS}
    case["verification_record_joint_trace_points"] = True
    remaining_times = affine_times(case, dt)
    flow = case["flow"]
    for material in flow["materials"]:
        material["volume_source_per_s"] = 0.
    head = dict(head_m=2., reference_point_m=[0., 0., 0.], gradient=[0., 0., slope])
    flow["initial_head"] = head
    flow["analysis"]["periods"][0]["boundaries"] = [
        dict(face_set=side, condition=dict(kind="head", value=head)) for side in ("bottom", "top")
    ]
    flow["settings"].update(initial_time_step_s=dt, maximum_time_step_s=dt)
    ground = controls.read(controls.EX / "ground-joint-water-mc.json")
    ground["initialization"].update(water_level_m=None)
    ground["initialization"]["geostatic"]["hydrostatic_face_sets"] = []
    folder = out / f"affine-{slope}-{dt}"
    run = solve(folder, solver, "solve-one-way", mesh, ground, case)
    records = []
    coords = np.array(mesh["nodes"])
    for frame in run["frames"]:
        if frame["phase"] != "transient":
            continue
        state = frame["state"]
        seconds = state["hydraulic"]["time_s"]
        if seconds not in remaining_times:
            raise ValueError("unexpected or repeated physical time")
        remaining_times.remove(seconds)
        traces = state["fluid"]["joint_traces"]
        keys = {(t["joint_id"], t["face_index"], t.get("point_index")) for t in traces}
        if len(traces) != len(expected_points) or keys != expected_points:
            raise ValueError("incomplete, duplicate, or unexpected JOINT traces")
        forces = np.array(state["mechanical"]["joint_water"]["face_nodal_force_a_kn"])
        if forces.shape != (len(faces), 3, 3):
            raise ValueError("incomplete or unexpected JOINT face nodal forces")
        for t in traces:
            a, b = t["side_a"], t["side_b"]
            if seconds != t["time_s"] or a["face"] != b["face"] or a["owner"] == b["owner"]:
                raise ValueError("unpaired time or continuous face")
            z = t["elevation_m"]
            exact_head = 2. + slope*z
            pressure = 10.*(exact_head-z)
            face = t["face_index"]
            xyz = coords[faces[face]]
            area = np.linalg.norm(np.cross(xyz[1]-xyz[0], xyz[2]-xyz[0]))/2
            cell_z = coords[mesh["tetrahedra"][a["owner"]], 2].mean()
            outward = math.copysign(1., z-cell_z)
            q = -1e-6*slope*area*outward
            expected_force = np.array([0., 0., -outward*pressure*area/3])
            records.append(dict(time_s=seconds, face=face, point_index=t["point_index"],
                                head_reference_m=exact_head,
                                side_a_head_error_m=a["total_head_m"]-exact_head,
                                side_b_head_error_m=b["total_head_m"]-exact_head,
                                pressure_error_kpa=t["applied_pressure_kpa"]-pressure,
                                expected_a_outflow_m3_s=q, a_flow_error_m3_s=a["outward_flow_m3_s"]-q,
                                b_flow_error_m3_s=b["outward_flow_m3_s"]+q,
                                max_nodal_force_error_kn=float(np.max(np.abs(forces[face]-expected_force)))))
    if remaining_times:
        raise ValueError("incomplete physical-time series")
    budget = {}
    controls.record("continuous", mesh, run, folder, budget)
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
               for slope in AFFINE_SLOPES for dt in AFFINE_STEPS]
    controls.write(args.output / "affine.json", results)
