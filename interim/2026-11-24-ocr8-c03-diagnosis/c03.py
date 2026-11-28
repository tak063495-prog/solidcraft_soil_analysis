"""Cross the incoming stress history with the last accepted strain partition."""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "2026-11-24-independent-local-references"
sys.path.insert(0, str(SOURCE))
from c03_point import EYE, LABELS, PointReturn, tensor


def voigt(s):
    return s[[0, 1, 2, 0, 1, 2], [0, 1, 2, 1, 2, 0]].tolist()


def trace(law, start, de, count):
    stress = start.copy()
    for _ in range(count):
        old = stress
        stress, _, label = law.step(old, de / count)
    last = de[0] / count
    trial = old[0] + law.g2 * last + law.lame * np.trace(last) * EYE
    eig = np.linalg.eigvalsh(trial)
    multiplier = (eig @ law.normal - law.k) / (law.normal @ law.shift)
    face = eig - multiplier * law.shift
    width = 1e-10 * max(abs(eig).max(), law.c, 1)
    return stress, {
        "label": LABELS[label.item()],
        "last_trial_principal_ascending_kpa": eig.tolist(),
        "trial_yield_kpa": float(eig @ law.normal - law.k),
        "face_candidate_upper_gap_kpa": float(face[2] - face[1]),
        "corrected_upper_gap_kpa": float(np.diff(np.linalg.eigvalsh(stress[0]))[-1]),
        "width_kpa": float(width),
        "signed_gap_over_width": float((face[2] - face[1]) / width),
    }


def investigate(probe):
    data = json.loads(gzip.decompress((SOURCE / "c03-inputs.json.gz").read_bytes()))["m0c1"]
    archived = json.loads(gzip.decompress((SOURCE / "c03-history.json.gz").read_bytes()))["m0c1"]
    summary = json.loads((SOURCE / "c03.json").read_bytes())["routes"]["m0c1"]
    law = PointReturn(data["material"])
    results, requests, checks = [], [], []
    for gp in (108, 185):
        j = data["points"].index(gp)
        event, = summary[j]["mismatches"]
        seed = tensor([data["initial_stress_kpa"][j]])
        parents = {"native": seed, "fine": seed}
        for i, row in enumerate(data["frames"]):
            de = tensor([row["strain"][j]], engineering=True)
            counts = {"native": row["fem_substeps"][j], "fine": 128}
            if (row["stage"], row["progress"]) == (event["stage"], event["progress"]):
                break
            for key, n in counts.items():
                parents[key], _, code = step_many(law, parents[key], de, n)
                expect = row["fem_label"][j] if key == "native" else LABELS[archived["reference_label"][i][j]]
                if LABELS[code.item()] != expect:
                    raise ValueError((gp, i, key, expect))
        else:
            raise ValueError("missing classified event")
        experiments = []
        for key, parent in parents.items():
            for which, count in counts.items():
                final, obs = trace(law, parent, de, count)
                if key == which and obs["label"] != event["fem" if key == "native" else "reference"]:
                    raise ValueError((gp, key, obs))
                experiments.append(dict(parent=key, increment_partition=which, count=count, **obs))
                # MC is perfectly plastic. The stress alone enters the predictor and return.
                old = dict(effective_stress_compression_positive_kpa=voigt(parent[0]),
                           void_ratio=1., plastic_strain_compression_positive=[0.]*6,
                           equivalent_plastic_strain=0., plastic_work_kpa=0., last_return="elastic")
                requests.append(dict(route="m0c1", increment=i, gp=gp, parameters=data["material"],
                                     old=old, strain=row["strain"][j], count=count))
                checks.append((experiments[-1], final[0]))
        results.append(dict(gp=gp, stage=row["stage"], progress=row["progress"],
                            prior_frames_replayed=i,
                            parent_stresses_kpa={k:voigt(v[0]) for k,v in parents.items()},
                            experiments=experiments))
    with tempfile.TemporaryDirectory() as folder:
        inp, out = (Path(folder) / name for name in ("input", "output"))
        inp.write_text(json.dumps(requests))
        subprocess.run([probe, inp, out], check=True)
        for (obs, final), native in zip(checks, json.loads(out.read_bytes()), strict=True):
            state = native["state"]
            obs["rust_label"] = state["last_return"]
            stress = tensor(state["effective_stress_compression_positive_kpa"])
            obs["rust_corrected_upper_gap_kpa"] = float(np.diff(np.linalg.eigvalsh(stress))[-1])
            obs["rust_stress_max_abs_delta_kpa"] = float(np.max(np.abs(stress - final)))
            if obs["rust_label"] != obs["label"]:
                raise ValueError(obs)
    return results


def step_many(law, stress, de, count):
    piece = de / count
    for _ in range(count):
        stress, plastic, label = law.step(stress, piece)
    return stress, plastic, label


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("probe", type=Path)
    args = p.parse_args()
    print(json.dumps(investigate(args.probe.resolve()), indent=2, allow_nan=False))
