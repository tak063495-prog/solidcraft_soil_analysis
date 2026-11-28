"""Slip, unload and transfer two independent Coulomb pad reactions to a shared-tip patch."""

import argparse
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
INPUT = HERE.parent / "2026-11-23-impl-root-causes-a-nc/a02-input"
sys.path.insert(0, str(HERE.parents[1] / "tools"))
from verify_joint_coupling import contact


def solve(probe, folder):
    folder.mkdir(parents=True, exist_ok=False)
    mesh, case = [json.loads((INPUT / (k + ".json")).read_bytes()) for k in ("mesh", "case")]
    (seam,) = case["joint_system"]["joints"]
    law = seam["material"]
    kn, kt, mu = (
        law[k] for k in ("normal_stiffness_kpa_m", "tangential_stiffness_kpa_m", "friction_coefficient")
    )
    system = SimpleNamespace(
        kn=np.array([kn]), kt=np.array([kt]), mu=np.array([mu]), normal=np.array([[0.0, 0.0, 1.0]])
    )
    xyz = np.array(mesh["nodes"])
    forces = np.zeros_like(xyz)
    pads = []
    facet_fields = {}
    # Independently displacement-controlled donor pads. Distinct donor DOFs
    # permit a jump in shear across the edge of the receiving facets.
    for face, pressure, direction, peak, unload in zip(
        seam["crack_tip"]["faces"],
        (60.0, 40.0),
        ((1.0, 0.0), (0.0, -1.0)),
        (0.00015, 0.00011),
        (0.000025, 0.00002),
        strict=True,
    ):
        tangent = np.array([*direction, 0.0])
        closed = np.array([0.0, 0.0, -pressure / kn])
        old_jump = old_force = np.zeros((1, 3))
        history = []
        for end in (closed, closed + peak * tangent, closed + (peak - unload) * tangent):
            tr, _, slip = contact(system, old_jump, old_force, end[None, :], 1)
            old_force, old_jump = tr[:, 0], end[None, :]
            history.append(
                dict(
                    jump_m=end.tolist(),
                    traction_kpa=old_force[0].tolist(),
                    plastic_slip_increment_m=float(slip.item()),
                )
            )
        expected = (mu * pressure - kt * unload) * tangent + [0.0, 0.0, -pressure]
        assert (
            history[1]["plastic_slip_increment_m"] > 0 and history[2]["plastic_slip_increment_m"] == 0
        )
        assert 0 < np.linalg.norm(expected[:2]) < mu * pressure
        # Analytic elastic unloading from the sliding cap supplies a second check.
        np.testing.assert_allclose(old_force[0], expected, rtol=0, atol=1e-10)
        facet_fields[face["owner_a"]] = expected
        a, b, c = xyz[face["face_a"]]
        area = np.linalg.norm(np.cross(b - a, c - a)) / 2
        for side, sign in (("a", -1), ("b", 1)):
            mids = mesh["tetrahedra10"][face["owner_" + side]][4:]
            for i, j in combinations(face["face_" + side], 2):
                (node,) = [n for n in mids if np.array_equal(xyz[n], (xyz[i] + xyz[j]) / 2)]
                forces[node] += sign * area * old_force[0] / 3
        pads.append(
            dict(
                owner_a=face["owner_a"],
                area_m2=area,
                history=history,
                reference_kpa=expected.tolist(),
                cone_margin_kpa=mu * pressure - np.linalg.norm(expected[:2]),
            )
        )
    assert np.any(forces) and len(facet_fields) == 2
    request, output = folder / "request.json", folder / "native.json"
    request.write_text(json.dumps(dict(operation="birth", residual_kn=forces.ravel().tolist())))
    subprocess.run(
        [probe.resolve(), INPUT / "mesh.json", INPUT / "case.json", request, output], check=True
    )
    run = json.loads(output.read_bytes())
    out = run["result"]
    reconstructed = np.asarray(out["nodal_force_kn"]).reshape(forces.shape)
    force_error = float(np.abs(reconstructed - forces).max())
    assert run["accepted"] and out["exact_repeat"] and force_error <= 1e-10
    errors, margins, points = [], [], []
    for p, state in zip(run["layout"], out["state"]["points"], strict=True):
        q = np.asarray(state["birth"]["initial_traction_kpa"])
        np.testing.assert_array_equal(q, state["traction_kpa"])
        n = np.asarray(p["normal"])
        normal = q @ n
        margin = -mu * normal - np.linalg.norm(q - normal * n)
        assert normal < 0 and margin > 0
        errors.append(float(np.linalg.norm(q - facet_fields[p["owner_a"]])))
        margins.append(float(margin))
        points.append(
            dict(
                face=p["face_index"],
                owner_a=p["owner_a"],
                owner_b=p["owner_b"],
                position_m=p["position_m"],
                reference_kpa=facet_fields[p["owner_a"]].tolist(),
                native_kpa=q.tolist(),
                cone_margin_kpa=margin,
            )
        )
    assert max(errors) < 1e-7  # existing vector-traction regression, tests/joint_law.rs
    return dict(
        accepted=True,
        nodes=len(forces),
        gp=len(errors),
        pads=pads,
        points=points,
        maximum_nodal_component_error_kn=force_error,
        maximum_gp_vector_error_kpa=max(errors),
        minimum_gp_cone_margin_kpa=min(margins),
        reference_force_kn=forces.tolist(),
        native_force_kn=reconstructed.tolist(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe", type=Path)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    result = solve(args.probe, args.folder)
    (args.folder / "comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
