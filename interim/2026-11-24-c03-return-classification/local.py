"""Replay physical inputs at the published disagreement points, without advancing the FEM path."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-23-ac-scope-decisions"))
from classifications import Model, interface, mc, native_update, read


def trace(mat, old, strain, count, native_probe=None, resolve=False):
    de = np.array(strain) / count
    state = old
    substeps = []
    for _ in range(count):
        if native_probe is None:
            state = mc.step(mat, state, de, resolve_edges=resolve, trace=substeps)
        else:
            state = native_update(native_probe, state, mat, de)
            sig = state["effective_stress_compression_positive_kpa"]
            substeps.append(
                dict(
                    label=state["last_return"],
                    principal_kpa=np.linalg.eigvalsh(mc.tensor(sig))[::-1].tolist(),
                )
            )
    return substeps


def run(matrix, probe):
    previous = read(HERE.parent / "2026-11-23-ac-scope-decisions/results/c03.json")
    results = []
    for route in ("m0c1", "m1c0", "m1c1"):
        model = Model(
            read(matrix / "mesh.json"),
            read(matrix / "ground.json"),
            read(matrix / (route + "-native/stages.json")),
        )
        points = {
            (r["gp"], r["route"].split("-")[-1]): r
            for r in previous["points"]
            if r["route"].startswith(route + "-")
        }
        for (gp, engine), own in points.items():
            other = points[gp, "native" if engine == "independent" else "independent"]
            result = dict(
                route=route,
                gp=gp,
                engine=engine,
                stage=own["stage"],
                progress=own["progress"],
                original_label=own["label"],
            )
            if own["field"] == "contact_status":
                n = model.system.normal[gp]
                jump = np.array(own["previous_jump_m"])
                force = np.array(own["previous_traction_kpa"])
                du = np.array(own["target_jump_m"]) - jump
                dq = np.array(other["previous_traction_kpa"]) - force
                # Exchange tangential traction and normal closure independently.
                for key, start, traction in [
                    ("own", jump, force),
                    ("other_parent_shear", jump, force + dq - n * (dq @ n)),
                    ("other_parent_closure", jump + n * ((other["previous_jump_m"] - jump) @ n), force),
                ]:
                    replay = interface(
                        model.system, gp, start, traction, start + du, own["count"], own["label"]
                    )
                    result[key + "_excess_kpa"] = replay["trial_minus_cap_kpa"]
            else:
                mat = model.mats[gp]
                old = own["previous_state"]
                de = own["total_compressive_strain"]
                count = own["count"]
                local_probe = probe if engine == "native" else None
                original = trace(mat, old, de, count, local_probe)
                assert original[-1]["label"] == own["label"]
                result["original_substeps"] = original
                result["resolved_substeps"] = (
                    original if local_probe else trace(mat, old, de, count, resolve=True)
                )
                alternate = trace(mat, old, de, count, None if local_probe else probe, resolve=True)
                result["same_input_other_solver_labels"] = [s["label"] for s in alternate]
                # MC moduli do not depend on accumulated plastic strain or void ratio.
                # Exchange just the six parent stress components, then just the total strain increment.
                for key, parent, strain, subdivisions in [
                    (
                        "other_parent_stress",
                        dict(
                            old,
                            effective_stress_compression_positive_kpa=other["previous_state"][
                                "effective_stress_compression_positive_kpa"
                            ],
                        ),
                        de,
                        count,
                    ),
                    ("other_strain", old, other["total_compressive_strain"], count),
                    ("other_count", old, de, other["count"]),
                ]:
                    history = trace(mat, parent, strain, subdivisions, local_probe, resolve=True)
                    result[key + "_labels"] = [s["label"] for s in history]
            results.append(result)
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("matrix", "probe", "output"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    with a.output.open("x", encoding="utf-8") as f:
        json.dump(run(a.matrix, a.probe.resolve()), f, indent=2, allow_nan=False)
        f.write("\n")
