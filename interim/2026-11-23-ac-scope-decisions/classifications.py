"""Point-local replay and strict-gate counterfactual at the published MC disagreements."""

import argparse
import hashlib
import heapq
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "2026-11-23-impl-root-causes-c"
sys.path.insert(0, str(PREVIOUS))
from c03 import Model, RETURNS, CONTACT, STAGES, measure, read
import c03diag_mc as mc
from verify_joint_coupling import contact

MATERIAL = "material_return"
CONTACT_STATE = "contact_status"


def stream(path):
    with path.open(encoding="utf-8") as f:
        for line in f:
            yield json.loads(line)


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def coordinate(row):
    return STAGES.index(row["stage"]), row["progress"]


def point_id(row):
    return (*coordinate(row), row["field"], row["gp"])


def principal_diagnostics(mat, previous, strain, returned):
    young, nu = mat["young_kpa"], mat["poisson"]
    shear = young / (2 * (1 + nu))
    stiffness = np.full((3, 3), young * nu / ((1 + nu) * (1 - 2 * nu))) + 2 * shear * np.eye(3)
    trial = np.asarray(previous) + np.r_[stiffness @ strain[:3], shear * strain[3:]]
    s = np.linalg.eigvalsh(mc.tensor(trial))[::-1]
    end = np.linalg.eigvalsh(mc.tensor(returned))[::-1]
    cohesion, phi, psi = mc.strength(mat)
    n = np.array([1 - np.sin(phi), 0, -1 - np.sin(phi)])
    flow = np.array([1 - np.sin(psi), 0, -1 - np.sin(psi)])
    excess = n @ s - 2 * cohesion * np.cos(phi)
    shift = stiffness @ flow
    face = s - shift * (max(excess, 0) / (n @ shift))
    tol = 1e-10 * max(max(abs(s)), mat["cohesion_kpa"], 1.0)
    return dict(
        trial_principal_kpa=s.tolist(),
        trial_yield_kpa=excess,
        native_elastic_threshold_kpa=tol,
        single_face_max_mid_kpa=face[0] - face[1],
        single_face_mid_min_kpa=face[1] - face[2],
        returned_principal_kpa=end.tolist(),
        returned_max_mid_kpa=end[0] - end[1],
        returned_mid_min_kpa=end[1] - end[2],
    )


def native_update(probe, old, mat, de):
    with TemporaryDirectory() as tmp:
        req, output = (Path(tmp) / s for s in ("request", "result"))
        req.write_text(
            json.dumps(
                [
                    dict(
                        route="point",
                        increment=0,
                        gp=0,
                        old=old,
                        parameters=mat,
                        strain=de.tolist(),
                        count=1,
                    )
                ]
            )
        )
        subprocess.run([probe, req, output], check=True)
        return read(output)[0]["state"]


def material(old, de, count, target, label, mat, probe):
    fraction = np.asarray(de) / count
    state = old
    for _ in range(count):
        previous = state["effective_stress_compression_positive_kpa"]
        state = native_update(probe, state, mat, fraction) if probe else mc.step(mat, state, fraction)
    end = state["effective_stress_compression_positive_kpa"]
    error = float(max(abs(np.array(end) - target)))
    assert state["last_return"] == label and error == 0.0, (label, state["last_return"], error)
    return dict(
        label=label,
        count=count,
        replay_stress_error_kpa=error,
        previous_state=old,
        total_compressive_strain=np.asarray(de).tolist(),
        **principal_diagnostics(mat, previous, fraction, end),
    )


def interface(system, i, old_target, old_q, end, count, label):
    sub = SimpleNamespace(**{k: getattr(system, k)[i : i + 1] for k in ("normal", "kn", "kt", "mu")})
    start = np.array(old_target)[None, :]
    end = np.array(end)[None, :]
    q = np.array(old_q)[None, :]
    if count > 1:
        before = start + (end - start) * (count - 1) / count
        q = contact(sub, start, q, before, count - 1)[0][:, 0]
    else:
        before = start
    n = sub.normal[0]
    assert before[0] @ n < 0 and end[0] @ n < 0
    delta = (end - before)[0]
    trial = q[0] - n * (q[0] @ n) + sub.kt[0] * (delta - n * (delta @ n))
    cap = -sub.mu[0] * sub.kn[0] * (end[0] @ n)
    return dict(
        label=label,
        count=count,
        gap_m=end[0] @ n,
        previous_jump_m=list(old_target),
        previous_traction_kpa=list(old_q),
        target_jump_m=end[0].tolist(),
        trial_shear_kpa=np.linalg.norm(trial),
        coulomb_limit_kpa=cap,
        trial_minus_cap_kpa=np.linalg.norm(trial) - cap,
    )


def material_state(f, gp):
    columns = {
        "stress": "effective_stress_compression_positive_kpa",
        "plastic_strain": "plastic_strain_compression_positive",
        "void": "void_ratio",
        "plastic_work": "plastic_work_kpa",
        "equivalent_plastic_strain": "equivalent_plastic_strain",
    }
    return {
        **{key: f[col][gp] for col, key in columns.items()},
        "last_return": RETURNS[f[MATERIAL][gp]],
    }


def collect_native(root, route, model, targets, probe):
    folder = root / (route + "-trace")
    assert sha(folder / "accepted.jsonl") == sha(root / (route + "-native") / "raw/accepted.jsonl")
    labels = {}
    observations = {}
    previous = None
    for curr, d in zip(
        stream(folder / "accepted.jsonl"), stream(folder / "constitutive.jsonl"), strict=True
    ):
        key = d["stage"], d["progress"]
        assert curr["stage_index"] == key[0] and key not in labels
        assert curr["stage_progress"] == key[1] or (key == (0, 1.0) and curr["stage_progress"] == 0)
        mats = curr["material_states"]
        joints = curr["joints"]["points"]
        labels[key] = {
            MATERIAL: [p["last_return"] for p in mats],
            CONTACT_STATE: [p["status"] for p in joints],
        }
        for field, i in targets.get(key, ()):
            if field == MATERIAL:
                p = mats[i]
                value = material(
                    previous["material_states"][i],
                    d["strain_increment"][i],
                    d["partitions"][i],
                    p["effective_stress_compression_positive_kpa"],
                    p["last_return"],
                    model.mats[i],
                    probe,
                )
            else:
                old = previous["joints"]["points"][i]
                value = interface(
                    model.system,
                    i,
                    old["jump_m"],
                    old["traction_kpa"],
                    joints[i]["jump_m"],
                    d["partitions"][len(mats) + i],
                    joints[i]["status"],
                )
            if field == CONTACT_STATE:
                assert value["label"] == ("slip" if value["trial_minus_cap_kpa"] > 0 else "stick")
            observations[(*key, field, i)] = dict(**value, relative_force=d["relative_force_residual"])
        previous = curr
    return labels, observations


def collect_independent(root, route, model, targets):
    folder = root / (route + "-independent")
    assert read(folder / "result.json")["complete"]
    rows = heapq.merge(
        stream(folder / "fields.jsonl"),
        stream(folder / "substeps.jsonl"),
        key=lambda v: coordinate(v["row"]),
    )
    labels = {}
    observations = {}
    previous = None
    for item in rows:
        row, f = item["row"], item["fields"]
        key = coordinate(row)
        if (32 * key[1]).is_integer():
            assert key not in labels
            labels[key] = {
                MATERIAL: [RETURNS[j] for j in f[MATERIAL]],
                CONTACT_STATE: [CONTACT[j] for j in f[CONTACT_STATE]],
            }
        for field, i in targets.get(key, ()):
            if field == MATERIAL:
                du = (np.array(f["u"]) - previous["u"]).ravel()
                de = -model.strains(du)[i]
                value = material(
                    material_state(previous, i),
                    de,
                    row["material_partitions"][i],
                    f["stress"][i],
                    labels[key][field][i],
                    model.mats[i],
                    None,
                )
            else:
                value = interface(
                    model.system,
                    i,
                    previous["jump"][i],
                    previous["traction"][i],
                    f["jump"][i],
                    row["contact_partitions"][i],
                    labels[key][field][i],
                )
            if field == CONTACT_STATE:
                assert value["label"] == ("slip" if value["trial_minus_cap_kpa"] > 0 else "stick")
            observations[(*key, field, i)] = dict(
                **value, relative_force=row["relative_force_residual"]
            )
        previous = f
    return labels, observations


def tight_labels(folder):
    run = read(folder / "run.json")
    labels = {}
    residuals = {}
    attempts = [a for a in run["geostatic"]["attempts"] + run["attempts"] if a["accepted"]]
    for state, a in zip(stream(folder / "accepted.jsonl"), attempts, strict=True):
        key = state["stage_index"], a["progress_to"]
        assert state["stage_progress"] == key[1] or (key == (0, 1.0) and state["stage_progress"] == 0)
        assert key not in labels and a["relative_force_residual"] <= 1e-10
        labels[key] = {
            MATERIAL: [p["last_return"] for p in state["material_states"]],
            CONTACT_STATE: [p["status"] for p in state["joints"]["points"]],
        }
        residuals[key] = a["relative_force_residual"]
    return labels, residuals


def differences(a, b):
    if a.keys() != b.keys():
        raise ValueError("different accepted coordinates")
    return {
        (*key, field, i): (left, right)
        for key, fields in a.items()
        for field, v in fields.items()
        for i, (left, right) in enumerate(zip(v, b[key][field], strict=True))
        if left != right
    }


def run(root, probe):
    comparisons = {p.stem: read(p) for p in sorted((PREVIOUS / "results").glob("c03-*.json"))}
    targets = {}
    for v in comparisons.values():
        for r in v["mismatches"]:
            targets.setdefault(coordinate(r), set()).add((r["field"], r["gp"]))
    groups = {}
    tight = {}
    gate_summary = {}
    samples = []
    results = {}
    inputs = {}
    cross = []
    full_fields = {}
    for route in ("m0c1", "m1c0", "m1c1"):
        model = Model(
            read(root / "mesh.json"),
            read(root / "ground.json"),
            read(root / (route + "-native") / "stages.json"),
        )
        t, rt = tight_labels(root / (route + "-tight"))
        for engine in ("native", "independent"):
            if engine == "native":
                labels, obs = collect_native(root, route, model, targets, probe)
            else:
                labels, obs = collect_independent(root, route, model, targets)
            name = route + "-" + engine
            assert len(labels) == 224 and labels.keys() == t.keys()
            groups[name] = labels
            inputs[name] = obs
            for (stage, progress, field, i), row in sorted(obs.items()):
                samples.append(
                    dict(
                        route=name,
                        stage=STAGES[stage],
                        progress=progress,
                        field=field,
                        gp=i,
                        **row,
                        tight_native_label=t[stage, progress][field][i],
                        tight_native_relative_force=rt[stage, progress],
                    )
                )
        base_run = read(root / (route + "-native") / "raw/run.json")
        base_force = [
            a["relative_force_residual"]
            for a in base_run["geostatic"]["attempts"] + base_run["attempts"]
            if a["accepted"]
        ]
        assert len(base_force) == len(rt) == 224
        gate_summary[route] = dict(
            original_above_one_e_minus_ten=sum(x > 1e-10 for x in base_force),
            original_max=max(base_force),
            tight_max=max(rt.values()),
            changed_native_labels=len(differences(t, groups[route + "-native"])),
        )
        tight[route] = [
            dict(stage=STAGES[s], progress=p, field=f, gp=i, native=v[0], independent=v[1])
            for (s, p, f, i), v in differences(t, groups[route + "-independent"]).items()
        ]
    for k, v in comparisons.items():
        diff = differences(groups[v["left"]], groups[v["right"]])
        assert diff == {point_id(r): (r["left"], r["right"]) for r in v["mismatches"]}, k
        results[k] = len(diff)
        measures = measure(
            root, root / v["left"], root / v["right"], root / "m1c1-native/raw/run.json", v["axis"]
        )
        assert measures["accepted"] == 224 and measures["full_fields"].keys() == v["full_fields"].keys()
        full_fields[k] = dict(
            accepted=224,
            fields=measures["full_fields"],
            identical_to_previous=measures["full_fields"] == v["full_fields"],
        )
        for key in diff:
            stage, progress, field, i = key
            for own, other in ((v["left"], v["right"]), (v["right"], v["left"])):
                a, b = inputs[own][key], inputs[other][key]
                record = dict(
                    pair=k,
                    stage=STAGES[stage],
                    progress=progress,
                    field=field,
                    gp=i,
                    parent_route=own,
                    other_route=other,
                    original_label=a["label"],
                    original_count=a["count"],
                    other_count=b["count"],
                )
                if field == MATERIAL:
                    mat = model.mats[i]
                    for name, count, use_native in (
                        ("other_law", a["count"], not own.endswith("native")),
                        ("other_partition", b["count"], own.endswith("native")),
                    ):
                        state = a["previous_state"]
                        de = np.array(a["total_compressive_strain"]) / count
                        for _ in range(count):
                            state = (
                                native_update(probe, state, mat, de)
                                if use_native
                                else mc.step(mat, state, de)
                            )
                        record[name] = state["last_return"]
                else:
                    # Hold this route's parent traction and jump fixed while
                    # imposing the displacement increment from the other route.
                    delta = np.array(b["target_jump_m"]) - b["previous_jump_m"]
                    end = np.array(a["previous_jump_m"]) + delta
                    value = interface(
                        model.system,
                        i,
                        a["previous_jump_m"],
                        a["previous_traction_kpa"],
                        end,
                        a["count"],
                        a["label"],
                    )
                    record["other_jump_increment_trial_minus_cap_kpa"] = value["trial_minus_cap_kpa"]
                cross.append(record)
    identities = {
        "mesh": sha(root / "mesh.json"),
        "ground": sha(root / "ground.json"),
        "native_probe": read(root / "identity.json")["solver_sha256"],
        "local_native_return_probe": sha(probe),
        "native_cases": {
            route: sha(root / (route + "-native/stages.json")) for route in gate_summary
        },
    }
    return dict(
        input_and_binary_sha256=identities,
        original_mismatch_counts=results,
        force_gates=gate_summary,
        tight_native_vs_independent=tight,
        local_exchanges=cross,
        full_fields=full_fields,
        points=samples,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("matrix", "return_probe", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    result = run(args.matrix, args.return_probe.resolve())
    with args.output.open("x") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
