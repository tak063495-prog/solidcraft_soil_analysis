"""Apply the approved T10 criteria to named probes and physical-time accounts."""

import argparse
import hashlib
import json
import math
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T8 = HERE.parent / "2026-11-22-impl-t8-acceptance-criteria"
T2 = HERE.parent / "2026-11-20-impl-t2-soil-birth-hydraulic-init"
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(HERE.parent / "2026-11-21-impl-t6-offset-kink-joint"))
from t1_water_balance import LIMITS, audit
from t6_patch import patch
from t6_reference import shape, solve as multiplier_solution
from t6_study import compare as compare_multipliers, run as run_joint


def relative(value, reference, width):
    value = float(value)
    if not math.isfinite(value) or not math.isfinite(reference) or reference == 0:
        raise ValueError(
            "relative comparison requires a finite nonzero reference and finite sample"
        )
    return dict(
        measured=value,
        reference=reference,
        relative_width=width,
        signed_relative=(value - reference) / abs(reference),
        passed=abs(value - reference) <= width * abs(reference),
    )


def sealed_measurement(name):
    expected = dict(
        line.split(maxsplit=1)[::-1]
        for line in (T8 / "observations.sha256").read_text().splitlines()
    )
    with zipfile.ZipFile(T8 / "observations.zip") as archive:
        data = archive.read(name)
    if hashlib.sha256(data).hexdigest() != expected[name]:
        raise ValueError(f"changed T8 measurement {name}")
    return json.loads(data)


def shells():
    rows = []
    for source in sealed_measurement("b01/standard/measurements.json"):
        value = source.get("displacement")
        row = {k: source[k] for k in ("case", "n", "distorted", "diagnostic")}
        row.update(measured=value, passed=None)
        if source["case"] == "pinched_cylinder":
            if value is not None:
                row.update(relative(value, -1.8248e-5, 0.03))
        elif source["case"] == "scordelis_lo":
            row["theory_references"] = {
                "conventional": -0.3024,
                "MR_signed_interval": [-0.30204, -0.30200],
                "KL_approximate": -0.3006,
            }
        else:
            row["theory_references"] = {
                "conventional": 0.0940,
                "Ko_table_15": 0.093,
                "KL": 0.093597,
            }
        rows.append(row)
    return rows


def cantilever(binary):
    output = json.loads(subprocess.check_output([binary]))
    layers = []
    for layer in output["layers"]:
        samples = layer["root_samples"]
        interpolation = np.array([shape(*p["rs"])[0] for p in samples])
        projections = []
        for sample in samples:
            # Rotate to the test's x axis, about the shared reference plane.
            axes_x = np.array(sample["mid_axes"])[:2, 0]
            mx, my, mxy = sample["m_attachment_kn"]
            moment = axes_x @ [[mx, mxy], [mxy, my]] @ axes_x
            projections.append([axes_x @ sample["q_kn_m"], moment])
        qx, mxx = 1000 * np.linalg.solve(interpolation, projections)[0]
        layers.append(
            dict(
                offset_m=layer["offset_m"],
                qx_n_m=relative(qx, -500.0, 0.001),
                mxx_nm_m=relative(mxx, 5000.0, 0.001),
            )
        )
    return dict(
        triangles=output["triangles"],
        tip_dz_m=relative(output["tip_dz_m"], -3.90625e-5, 0.005),
        layers=layers,
    )


def joints(binary, work):
    out = {}
    for law in ("rigid", "hinge", "spring"):
        case = patch(1, law)
        run = run_joint(binary, case, work / law)
        state = run["frames"][-1]["active"]["joint"]
        left, right = np.array_split(state["stations"], 2)
        angles = np.array(
            [
                np.dot(
                    a["line_axis"], np.subtract(b["rotation_rad"], a["rotation_rad"])
                )
                for a, b in zip(left, right, strict=True)
            ]
        )
        moments = np.array(
            [np.dot(b["line_axis"], b["moment_at_joint_kn_m"]) for b in right]
        )
        if law == "rigid":
            result = dict(
                max_angle_rad=float(max(abs(angles))),
                passed=bool((abs(angles) <= 1e-14).all()),
            )
        elif law == "hinge":
            result = dict(
                max_moment_kn_m=float(max(abs(moments))),
                passed=bool((abs(moments) <= 1e-10).all()),
            )
        else:
            mass = np.array([[4, 2, -1], [2, 16, 2], [-1, 2, 4]]) / 30
            expected = -25 * mass @ angles
            energy = float(25 * angles @ mass @ angles / 2)
            actual = state["spring_energy_kn_m"]
            result = dict(
                max_moment_error_kn_m=float(max(abs(moments - expected))),
                energy_kn_m=actual,
                reference_energy_kn_m=energy,
                passed=bool(
                    np.allclose(moments, expected, atol=1e-10, rtol=1e-7)
                    and round(abs(actual - energy), 13) == 0
                ),
            )
        global_comparison = compare_multipliers(case, run)
        ref = multiplier_solution(case)
        global_comparison["independent_compliance_kn_m"] = float(ref["load"] @ ref["q"])
        for cut in global_comparison["cuts"].values():
            e = next(
                f["offset_m"]
                for f in case["shell"]["faces"]
                if f["group"] == cut["group"]
            )
            cut["M_attachment_difference"] = (
                np.array(cut["M_difference"]) - e * np.array(cut["N_difference"])
            ).tolist()
        result["multipliers"] = global_comparison
        out[law] = result
    return out


def edf_force_points():
    # ssnp162d.comm lines 29-60 and 456-488. EDF regressions are reported
    # separately from this solver's cohesionless patch (the EDF law includes adhesion).
    adhesion, mu, kn, young = 1e5, 0.5, 1e12, 3e10
    strength = adhesion / mu
    closure = (strength / kn + strength / young) / 3
    force = kn * closure * mu + adhesion
    return [
        dict(inst=t, **relative(calc, exact, 0.001))
        for t, calc, exact in [
            (5, 100041.099959, adhesion),
            (15, -1244566.5998, -force),
            (20, 1244646.51058, force),
        ]
    ]


def published_force_points():
    """Compare the EDF regression values to the printed model D force table."""
    table = (1.0e5, -1.2445666e6, 1.244647e6)
    return [
        dict(inst=row["inst"], **relative(row["measured"], printed, 0.001))
        for row, printed in zip(edf_force_points(), table, strict=True)
    ]


def contact_patch(binary):
    material = dict(
        normal_stiffness_kpa_m=1e6,
        tangential_stiffness_kpa_m=2e5,
        friction_coefficient=0.3,
    )
    pressure = 1000.0
    radius = material["friction_coefficient"] * pressure

    def path(initial_pressure, increments):
        case = dict(
            material=material,
            normal=[0, 0, 1],
            reference=[0, 0, 0],
            traction=[0, 0, -initial_pressure],
            increments=increments,
            maximum_increment_m=1e-5,
        )
        return json.loads(
            subprocess.check_output([binary], input=json.dumps(case).encode())
        )["states"]

    # Start at a physically admissible compression and load monotonically in x.
    positions = [0.001, 0.002, 0.01]
    states = path(pressure, [[x, 0, 0] for x in positions])[1:]
    rows = []
    for x, state in zip(positions, states, strict=True):
        elastic = radius / material["tangential_stiffness_kpa_m"]
        slip = max(0.0, x - elastic)
        expected_force = min(radius, x * material["tangential_stiffness_kpa_m"])
        work = radius * slip
        measured = state["plastic_dissipation_kpa_m"]
        rows.append(
            dict(
                x_m=x,
                state=state["status"],
                pressure_kpa=-state["traction_kpa"][2],
                force=relative(state["traction_kpa"][0], expected_force, 0.001),
                slip_m=state["accumulated_plastic_slip_m"],
                reference_slip_m=slip,
                work_kpa_m=measured,
                reference_work_kpa_m=work,
                work_passed=abs(measured - work) <= 0.001 * work,
            )
        )
    opening = path(0.0, [[0, 0, 0.002], [0.001, 0, 0.002]])[-1]
    return dict(
        points=rows,
        open_work_kpa_m=opening["plastic_dissipation_kpa_m"],
        open_state=opening["status"],
        open_passed=opening["status"] == "open"
        and opening["plastic_dissipation_kpa_m"] == 0.0,
    )


def water():
    placements = []
    for row in sealed_measurement("t8-water/placement.json"):
        level = row["level_m"]
        lower, upper = level - 0.5, 1.0 - level
        porosity = 0.7 / 1.7
        expected = dict(
            lower_volume_m3=lower,
            upper_volume_m3=upper,
            lower_water_m3=porosity * lower,
            upper_water_m3=porosity * 0.25 * upper,
            lower_weight_kn=20.0 * lower,
            upper_weight_kn=((1 - porosity) * 27 + porosity * 0.25 * 10) * upper,
        )
        errors = {key: row["measured"][key] - v for key, v in expected.items()}
        placements.append(
            dict(
                level_m=level,
                expected=expected,
                measured=row["measured"],
                difference=errors,
                passed=all(abs(x) <= 1e-12 for x in errors.values()),
            )
        )
    raw = (T2 / "evidence.zip").read_bytes()
    if (
        hashlib.sha256(raw).hexdigest()
        != json.loads((T2 / "acceptance.json").read_bytes())["evidence_zip_sha256"]
    ):
        raise ValueError("changed T2 hydraulic evidence")
    with zipfile.ZipFile(T2 / "evidence.zip") as archive:
        mesh, report = [
            json.loads(archive.read("controls/lifecycle/" + name + ".json"))
            for name in ("mesh", "result")
        ]
    (removal,) = report["excavations"]
    removed = removal["removed_hydraulic_points"]
    expected_cells = {int(i) for i in removal["removed_generations"]}
    expected_points = {(i, gp) for i in expected_cells for gp in range(27)}
    keys = {(p["element_index"], p["local_index"]) for p in removed}
    if len(expected_cells) != 6 or keys != expected_points or len(keys) != len(removed):
        raise ValueError("incomplete excavation quadrature")
    volume = math.fsum(p["weight_m3"] for p in removed)
    hand_volume = 1.0 * 1.0 * (1.0 - 0.5)
    excavation = dict(
        time_s=removal["time_s"],
        measured_m3=volume,
        hand_m3=hand_volume,
        absolute_error_m3=abs(volume - hand_volume),
        passed=abs(volume - hand_volume) <= 1e-12,
    )
    lifecycle = budget_profiles(mesh, report, audit)
    old = HERE.parent / "2026-11-22-wet-joint-transient-water"
    t1 = HERE.parent / "2026-11-20-impl-t1-wet-joint-richards" / "study.py"
    spec = importlib.util.spec_from_file_location("t10_t1_controls", t1)
    controls = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(controls)
    sys.path.insert(0, str(T8))
    from sealed import archive as open_parts

    with open_parts(old, "archive-sha256.json") as z:
        soil = budget_profiles(
            json.loads(z.read("data/final/hydraulic-mesh.json")),
            json.loads(z.read("data/final/hydraulic-h16/native.json")),
            controls.fixed_domain,
        )
    return dict(
        excavation=excavation,
        placements=placements,
        lifecycle=lifecycle,
        fixed_soil=soil,
        steady_domain_balance=None,
    )


def budget_profiles(mesh, report, auditor):
    profiles = {}
    for profile in LIMITS:
        balance = auditor(mesh, report, profile)
        accounts = [r for r in balance["rows"] if "passed" in r]
        global_rows = [r for r in accounts if r["kind"] == "global"]
        if not global_rows or balance["passed"] is None:
            raise ValueError("missing lifecycle inventory or physical-time balance")
        profiles[profile] = dict(
            passed=balance["passed"],
            account_rows=len(accounts),
            failures=sum(not r["passed"] for r in accounts),
            global_rows=global_rows,
        )
    return profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    suffix = ".exe" if os.name == "nt" else ""
    release = ROOT / "target/release"
    from t10_history import wet_history

    with tempfile.TemporaryDirectory() as tmp:
        report = {
            "B-01": shells(),
            "B-02": cantilever(
                release / "examples" / ("t10_layered_cantilever" + suffix)
            ),
            "B-03": joints(release / ("soil-fem-3d" + suffix), Path(tmp)),
            "C-04": {
                "edf_input_force_checks": edf_force_points(),
                "patch": contact_patch(
                    release / "examples" / ("joint_birth_points" + suffix)
                ),
                "history": wet_history(),
            },
            "A-03": water(),
        }
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    main()
