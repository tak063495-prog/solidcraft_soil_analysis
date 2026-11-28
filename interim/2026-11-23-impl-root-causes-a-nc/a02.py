"""Verify the compressive column transfer against a closed-form Tri6 load vector."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess

import numpy as np


INPUT = Path(__file__).resolve().with_name("a02-input")
FORCE_TOLERANCE_KN = 1e-10


def exact_load(mesh, ground, case):
    xyz = np.array(mesh["nodes"])
    pairs = case["joint_system"]["joints"][0]["crack_tip"]["faces"]
    seam_z = xyz[pairs[0]["face_a"][0], 2]
    upper_height = xyz[:, 2].max() - seam_z
    height = np.ptp(xyz[:, 2])
    top = next(b for b in case["boundary"]["constraints"]
               if b["node_set"] == "top" and b["direction"] == [0, 0, 1])
    delta = top["value_m"] * ground["initialization"]["geostatic"]["initial_boundary_factor"]
    mat = ground["materials"][0]
    # Lower half equilibrated before placement. Integrate the upper half's
    # triangular self-weight stress and impose total column shortening.
    normal = (mat["mohr_coulomb"]["young_kpa"] * delta
              - mat["unit_weight_kn_m3"] * upper_height**2 / 2) / height
    expected = np.zeros_like(xyz)
    for pair in pairs:
        for side, sign in (("a", -1), ("b", 1)):
            vertices = pair["face_" + side]
            assert all(xyz[i, 2] == seam_z for i in vertices)
            a, b, c = xyz[vertices]
            area = np.linalg.norm(np.cross(b-a, c-a)) / 2
            mids = mesh["tetrahedra10"][pair["owner_" + side]][4:]
            for i, j in itertools.combinations(vertices, 2):
                node, = [n for n in mids if np.array_equal(xyz[n], (xyz[i]+xyz[j])/2)]
                # Exact area integrals on affine Tri6: corners 0, midsides A/3.
                # Use physical IDs, including the nodes shared across the tip.
                expected[node, 2] += sign * normal * area / 3
    return normal, expected


def verify(run, mesh, ground, case):
    normal, expected = exact_load(mesh, ground, case)
    assert normal < 0 and np.any(expected)
    definition, = case["joint_system"]["joints"]
    mu = definition["material"]["friction_coefficient"]
    layout = run["joint_integration_points"]
    frames = [f["state"] for f in run["frames"] if f["name"] == "birth"]
    assert run["complete"] and frames and frames[-1]["stage_progress"] == 1
    assert len(layout) == 32
    rows = []
    for frame in frames:
        joints = frame["joints"]
        assert all(r["active"] and r["generation"] == 1 for r in joints["lifecycle"])
        for field in ("initial_traction_kpa", "traction_kpa"):
            force = np.zeros_like(expected)
            gp_error = 0.
            margin = np.inf
            for p, state in zip(layout, joints["points"], strict=True):
                q = np.array((state["birth"] if field == "initial_traction_kpa" else state)[field])
                n = np.array(p["normal"])
                qn = q @ n
                cone_margin = -mu * qn - np.linalg.norm(q - qn*n)
                assert qn <= 0 and cone_margin >= 0
                margin = min(margin, cone_margin)
                gp_error = max(gp_error, np.max(np.abs(q - [0, 0, normal])))
                f = np.outer(p["shape"], q) * p["weight_m2"]
                force[p["nodes_a"]] -= f
                force[p["nodes_b"]] += f
            error = float(np.max(np.abs(force - expected)))
            assert error <= FORCE_TOLERANCE_KN, (field, error)
            rows.append(dict(progress=frame["stage_progress"], field=field,
                             maximum_component_error_kn=error,
                             maximum_gp_error_kpa=gp_error, minimum_cone_margin_kpa=margin))
    return dict(limited_pass=True, normal_reference_kpa=normal, nodes=len(expected),
                components=expected.size, gp=len(layout), accepted_birth=sum(a["accepted"] and a["stage"] == "birth"
                                   for a in run["attempts"]), saved_birth_frames=len(frames), rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cli", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    files = [INPUT / (k + ".json") for k in ("mesh", "ground", "case")]
    run = args.output / "run.json"
    subprocess.run([args.cli.resolve(), "solve-stages", *files, run], check=True)
    result = verify(json.loads(run.read_bytes()), *(json.loads(p.read_bytes()) for p in files))
    (args.output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                                              encoding="utf-8")


if __name__ == "__main__":
    main()
