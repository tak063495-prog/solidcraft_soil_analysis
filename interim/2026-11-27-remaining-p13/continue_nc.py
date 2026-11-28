"""Continue each p12 checkpoint, recording the p10 prediction before every parent."""
import importlib.util
import gzip
import json
from pathlib import Path
import sys

import numpy as np
from plans import HERE, P12, save

sys.path.insert(0, str(P12))
spec=importlib.util.spec_from_file_location("p12_continuation",P12/"nc.py")
p12=importlib.util.module_from_spec(spec)
spec.loader.exec_module(p12)
advance,digest,old,predict=p12.advance,p12.digest,p12.old,p12.predict


def run(work, out):
    arc = json.loads((HERE / "design.json").read_text())["nc"]["common_arc_m"]
    prior = json.loads((P12 / "nc/comparison.json").read_text())
    count = prior["cumulative_including_replay"]
    probe = work / "probe/target/release/examples"
    model = old.nc_extended.Extended(work / "inputs", probe / "p5_homogeneous")
    assert len(model.free) == 348 and len(model.b) == 48
    paths = {}
    with (out / "global.jsonl").open("x", buffering=1) as log:
        def record(**kw):
            print(json.dumps(kw, allow_nan=False), file=log)

        for path in ("ordinary", "predicted"):
            data = old.read(P12 / "nc" / (path + "-end.json.gz"))
            parent = (data["active_material_states"], data["load"], np.array(data["secant"]),
                      np.array(data["displacement_m"]), data["next_width"])
            start = digest(parent)
            assert start == data["state_hash"] == prior["paths"][path]["end_hash"]
            used, accepted, length, chord = count, 0, 0., 0.
            stop = None
            movement = []
            with gzip.open(out / (path + "-local.jsonl.gz"), "xt") as local:
                while length < arc:
                    state, lam, direction, _, width = parent
                    key = digest(parent)
                    row = predict(model, probe, state, lam, direction, width, accepted + 1, local)
                    record(path=path, parent_hash=key, arc_before=length, frozen_prediction=row)
                    movement.append(dict(arc=length, load=lam, boundaries=row["partition_boundaries"],
                        roots=[[r["point"], r["left"]["partition"], r["left"]["width"], r["right"]["width"]]
                               for r in row["sign_brackets"]],
                        unresolved=[[r["point"], r["left"]["width"], r["right"]["width"]]
                                    for r in row["unresolved_brackets"]]))
                    if not row["unique_width"]:
                        stop = "p10 minimum unresolved"
                        break
                    w = min(width if path == "ordinary" else row["predicted"], arc - length)
                    result = None
                    while result is None:
                        if count >= min(2000, (2065 * 96 + 11) // 12):
                            stop = "attempt limit"
                            break
                        count += 1
                        result, it, why = old.nc_branch.correct(model, state, lam, direction, w, 1e-12)
                        assert digest(parent) == key
                        record(path=path, parent_hash=key, load=lam, width=w, accepted=result is not None,
                               cumulative=count, iterations=it, reason=why)
                        if result is None:
                            w /= 2
                            if w < model.settings["minimum_increment"]:
                                stop = "minimum increment"
                                break
                    if stop:
                        break
                    parent = advance(model, parent, w, result, it)
                    length += w
                    chord += float(np.linalg.norm(result[0]))
                    accepted += 1
                    record(path=path, accepted_parent=accepted, arc=length, chord=chord,
                           state_hash=digest(parent), load=parent[1], next_width=parent[-1])
            state, lam, direction, u, width = parent
            with gzip.open(out / (path + "-end.json.gz"), "xt") as f:
                json.dump(dict(active_material_states=state, load=lam, secant=direction.tolist(),
                               displacement_m=u.tolist(), next_width=width, state_hash=digest(parent)), f, allow_nan=False)
            paths[path] = dict(start_hash=start, end_hash=digest(parent), start_load=data["load"],
                               end_load=lam, delta_lambda=lam-data["load"], secant_lambda=direction[-1],
                               arc=length, chord=chord, accepted=accepted, rejected=count-used-accepted,
                               attempts=count-used, movement=movement, stop=stop)
            save(out / "comparison.json", dict(common_arc=arc, paths=paths, cumulative=count))


if __name__ == "__main__":
    run(Path(sys.argv[1]).resolve(), Path(sys.argv[2]))
