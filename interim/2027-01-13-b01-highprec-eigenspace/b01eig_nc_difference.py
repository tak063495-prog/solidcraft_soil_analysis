"""Measure the saved NC force split using only verified snapshot inputs.

Run with python -I. The official pass/fail test is run separately at both
commits. This read-only observer emits the signed component differences in kN.
"""

import argparse
from contextlib import closing
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
NC_REL = Path("interim/2026-11-10-nc-mean-stress-release")
ARCHIVE_SHA = "87434944b291e5958f8a98e3d36ff485a64adbe95be67f42d3fbbe2cee25f6f9"
NATIVE_SEED = "prepared/origin/prior/fresh-head-midpoint/job/recovery-checkpoint-0.json"
SEEDS = {NATIVE_SEED, "prepared/reference400/unload.json", "prepared/reference800/unload.json"}


def read(p):
    return json.loads(p.read_bytes())


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def extract_inputs(dest, manifest, folder):
    replay = load_module("b01eig_archived_reader", folder / "remeasure.py")
    archive_hash = hashlib.sha256()
    parts = []
    for entry in manifest["parts"]:
        part = folder / entry["name"]
        raw = part.read_bytes()
        if len(raw) != entry["size"] or hashlib.sha256(raw).hexdigest() != entry["sha256"]:
            raise ValueError(f"part identity differs: {part.name}")
        archive_hash.update(raw)
        parts.append(part)
    if archive_hash.hexdigest() != ARCHIVE_SHA or manifest["archive_sha256"] != ARCHIVE_SHA:
        raise ValueError("NC archive identity differs")
    # measure() discovers direct children with controls AND either a report or a watchdog.
    entries = manifest["entries"]
    cases = sorted(
        name for path in entries
        if path.count("/") == 1 and path.endswith("/controls.json")
        for name in [path.split("/")[0]]
        if f"{name}/report.json" in entries or f"{name}/external-watchdog.json" in entries
    )
    wanted = SEEDS | {c + "/controls.json" for c in cases}
    # make(mcc=True) reads the saved mesh, ground and stage definition.
    wanted.update({
        "sources/examples/joint-column-tet10-mesh.json",
        "sources/examples/ground-joint-mcc.json",
        "sources/examples/stages-joint-explicit-normal-column.json",
    })
    wanted.update(p for p in entries if p.startswith("sources/") and p.endswith(".py"))
    seen = set()
    with closing(replay.PartsReader(parts)) as reader, tarfile.open(fileobj=reader, mode="r|gz") as archive:
        for ent in archive:
            if ent.name not in wanted:
                continue
            path = Path(ent.name)
            if path.is_absolute() or ".." in path.parts or not ent.isfile() or ent.name in seen:
                raise ValueError(ent.name)
            raw = archive.extractfile(ent).read()
            expected = entries[ent.name]
            if len(raw) != expected["size"] or hashlib.sha256(raw).hexdigest() != expected["sha256"]:
                raise ValueError(f"member identity differs: {ent.name}")
            seen.add(ent.name)
            out = dest / path
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(raw)
    if not cases or seen != wanted:
        raise ValueError(f"missing diagnostic inputs: {sorted(wanted - seen)}")
    return cases


def compare(case, dest, measure, mesh):
    c = read(dest / case / "controls.json")
    active = np.array(c["active_elements"])
    if c["binding"]["map"] == "native":
        seed = read(dest / NATIVE_SEED)
        state = next(f["state"] for f in seed["run"]["frames"] if f["name"] == "unload")
        stress = np.array([state["material_states"][i]["effective_stress_compression_positive_kpa"]
                           for i in c["active_gp"]])
    else:
        seed = "reference800" if case.endswith("-pre800") else "reference400"
        stress = np.array(read(dest / "prepared" / seed / "unload.json")["material_state"])[active].reshape(-1, 14)[:, :6]
    if not np.array_equal(stress, c["split_source_stress_kpa"]):
        raise ValueError(f"pre-cut stress differs: {case}")
    bounded = c["binding"].get("force_decomposition") == "sign_preserving_bounded_mean"
    partition = measure.bounded_partition if bounded else measure.partition
    order = c["binding"]["order"]
    recomputed = partition(mesh, active, stress, np.array(c["holding_kn"]), order)
    differences = []
    keys = ("horizontal_holding_kn", "vertical_holding_kn")
    splits = ("mean", "complement") if order == "mean-first" else ("complement", "mean")
    for saved_key, split, new in zip(keys, splits, recomputed, strict=True):
        old = np.array(c[saved_key])
        if new.shape != old.shape or old.shape != (3 * len(mesh["nodes"]),):
            raise ValueError(f"holding force layout differs: {case}")
        for i in np.flatnonzero(new != old):
            differences.append({
                "saved_key": saved_key,
                "split": split,
                "dof": int(i),
                "node_zero_based": int(i // 3),
                "component": "xyz"[i % 3],
                "saved_kn": float(old[i]),
                "recomputed_kn": float(new[i]),
                "signed_difference_kn": float(new[i] - old[i]),
                "saved_hex": float(old[i]).hex(),
                "recomputed_hex": float(new[i]).hex(),
            })
    return {
        "case": case, "map": c["binding"]["map"], "order": order, "bounded": bounded,
        "differences": differences,
        "max_abs_difference_kn": max((abs(r["signed_difference_kn"]) for r in differences), default=0.0),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=ROOT, help="checkout supplying the sealed archive")
    args = parser.parse_args()
    if not sys.flags.isolated:
        parser.error("use python -I so the archived imports are isolated")
    folder = args.root.resolve() / NC_REL
    manifest = read(folder / "archive-sha256.json")
    # Loading the streaming reader imports live tools. Start the measurement
    # in a clean interpreter so every NC helper resolves to the sealed snapshot.
    with tempfile.TemporaryDirectory() as tmp:
        dest = Path(tmp)
        cases = extract_inputs(dest, manifest, folder)
        (dest / "case-names.json").write_text(json.dumps(cases), encoding="utf-8")
        subprocess.run(
            [sys.executable, "-I", "-X", "utf8", str(Path(__file__).resolve()), "--snapshot",
             str(dest), str(args.output.resolve())],
            check=True, cwd=dest,
        )


if __name__ == "__main__":
    if sys.argv[1:2] == ["--snapshot"]:
        dest, output = map(Path, sys.argv[2:])
        snapshot = dest / "sources"
        measure = load_module("b01eig_sealed_measure", snapshot / NC_REL / "measure.py")
        mesh, _, _ = measure.make(mcc=True, pc=100., preload_kpa=5., force=(1., .5, 0.))
        cases = [compare(case, dest, measure, mesh) for case in read(dest / "case-names.json")]
        report = {
            "archive_sha256": ARCHIVE_SHA,
            "numpy_version": np.__version__,
            "first_mismatching_case": next((r["case"] for r in cases if r["differences"]), None),
            "max_abs_difference_kn": max(r["max_abs_difference_kn"] for r in cases),
            "cases": cases,
        }
        with output.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        main()
