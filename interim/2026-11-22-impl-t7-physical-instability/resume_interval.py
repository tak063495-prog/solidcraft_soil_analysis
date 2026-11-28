"""Continue an unfinished diagnostic interval from its last accepted state."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import ijson
import numpy as np
from continuation import Controls, assess_path, identity, serial, spectrum, trace
from parallel import Model, Native, U, independent, save


def read_field(path, key, multiple=False):
    with gzip.open(path) as stream:
        items = ijson.items(stream, key, use_float=True)
        if not multiple:
            return next(items)
        last = None
        for last in items:
            pass
        if last is None:
            raise ValueError("no accepted checkpoint")
        return last


def resume(folder, label, bridge, source, source_sha, output):
    with source.open("rb") as f:
        digest = hashlib.file_digest(f, "sha256").hexdigest()
    if digest != source_sha:
        raise ValueError("archived trace digest differs")
    data = [
        json.loads((folder / (k + ".json")).read_bytes())
        for k in ("mesh", "ground", "case")
    ]
    model = Model(*data)
    meta = read_field(source, "environment")
    if read_field(source, "input_sha") != {
        k: identity(v) for k, v in zip(("mesh", "ground", "case"), data)
    }:
        raise ValueError("mesh, material, or stage input changed")
    if hashlib.sha256(bridge.read_bytes()).hexdigest() != meta["bridge_sha256"]:
        raise ValueError("evaluation binary differs")
    for name, expected in meta["sources"].items():
        if (
            hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            != expected
        ):
            raise ValueError("continuation source differs: " + name)
    if float(label.split("-")[1]) != read_field(source, "trials.item")["ds"] or int(
        label.split("-")[2]
    ) != read_field(source, "subdivision"):
        raise ValueError("arc or contact step differs")
    last = read_field(source, "last_point")
    restart = read_field(source, "restart")
    parent = last["state"]
    z = np.asarray(last["z"])
    engine, step, parts = label.split("-")
    if engine == "native":
        op = Native(
            bridge,
            folder,
            model,
            json.loads((folder / "run.json").read_bytes()),
            int(parts),
        )
        tol = model.case["settings"]["relative_force_tolerance"]
    else:
        with gzip.open(folder / "independent-last.json.gz") as f:
            record = json.load(f)
        path = record["stage_path"]
        path["holding"] = np.asarray(path["holding"])
        op = independent(model, path, int(parts))
        tol = 1e-10
        parent["u"] = np.asarray(parent["u"])
        parent["joint"] = {k: np.asarray(v) for k, v in parent["joint"].items()}
    try:
        # Bound storage per segment. Emit the last full checkpoint before
        # advancing so a process interruption cannot erase its parent.
        result = trace(
            op,
            parent,
            z,
            restart["payload"]["previous_direction"],
            Controls(float(step), 64, tol, 1e-10, 18, event_bisections=4),
            restart=restart,
            stop=lambda pts: pts[-1]["z"][-1] >= 0.125,
        )
        result["fold_candidates"] = assess_path(op, result["accepted"], tol)
        for record in result["critical"].values():
            if "iterations" not in record:
                continue
            ev = next(
                (
                    i["evaluation"]
                    for i in reversed(record["iterations"])
                    if "evaluation" in i
                ),
                None,
            )
            if ev:
                ev["spectrum"] = spectrum(ev["k"] / U, ev["q"], ev["residual"])
    finally:
        if engine == "native":
            op.close()
    result.update(
        parent_trace_sha256=source_sha,
        parent_s=last["s"],
        stage=read_field(source, "stage"),
        input_sha=read_field(source, "input_sha"),
        environment=meta,
        subdivision=int(parts),
    )
    save(output, result)
    return dict(
        stop=result["termination_reason"],
        start=z[-1],
        end=result["accepted"][-1]["z"][-1],
        accepted=len(result["accepted"]) - 1,
        attempts=len(result["trials"]),
        last=result["last_accepted"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    parser.add_argument("label")
    parser.add_argument("bridge", type=Path)
    parser.add_argument("source", type=Path)
    parser.add_argument("source_sha")
    parser.add_argument("output", type=Path)
    a = parser.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    print(
        json.dumps(
            resume(
                a.folder, a.label, a.bridge.resolve(), a.source, a.source_sha, a.output
            ),
            default=serial,
        )
    )
