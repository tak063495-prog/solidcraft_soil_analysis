"""Compare only matching accepted open coordinates and verify both real restarts."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

from branch import Model, Native, read, sha, unpack
from a01follow_binding import fields


def segments(folder, engine):
    direct, first, second = [unpack(folder / (engine + suffix + ".gz"))
                             for suffix in ("", "-4", "-resumed")]
    assert first["stop"] == "segment end" and second["parent_sha256"] == sha(folder / (engine + "-4.gz"))
    assert first["origins"] == second["origins"] == direct["origins"] == sha(folder / "origins.gz")
    assert first["points"] + second["points"] == direct["points"]
    assert first["trials"] + second["trials"] == direct["trials"]
    for key in ("last_accepted", "last_progress", "stop", "restart", "steps_used"):
        assert direct[key] == second[key]
    return direct


def compare(folder, bridge):
    native, reference = [segments(folder, engine) for engine in ("native", "independent")]
    model = Model(*(read(folder / (k + ".json")) for k in ("mesh", "ground", "case")))
    data = read(folder / "native.json")
    data["last_committed"] = unpack(folder / "origins.gz")["native"]
    op = Native(bridge.resolve(), folder, model, data, 1)
    common = []
    independent_by_load = {(p["z"][-1], p["factor"]): p for p in reference["points"]}
    parent = data["last_committed"]
    try:
        for n in native["points"]:
            state = n["state"]
            ref = independent_by_load.get((n["z"][-1], n["factor"]))
            if ref is not None:
                assert n["event"] == ref["event"]
                # Re-evaluate the same candidate against the actual accepted parent
                # to obtain physical reactions in the native ordering.
                ev = op(parent, np.asarray(n["z"]))
                target = np.asarray(state["applied_nodal_load_kn"]).ravel()
                reaction = np.asarray(ev.details["reaction"]).ravel()
                ordered = copy.deepcopy(state)
                ordered["joints"]["points"] = [state["joints"]["points"][i] for i in op.order]
                left = fields.native_fields(model, ordered, dict(internal=target+reaction, target=target), {})
                right = ref["state"]["fields"]
                row = {}
                for name, value in right.items():
                    x, y = left[name], np.asarray(value)
                    assert x.shape == y.shape and np.isfinite(x).all() and np.isfinite(y).all()
                    if name in fields.CATEGORIES:
                        row[name] = dict(mismatches=int(np.count_nonzero(x != y)))
                    else:
                        row[name] = dict(absolute_l2=float(np.linalg.norm(x-y)))
                assert row["generation"]["mismatches"] == row["contact_status"]["mismatches"] == 0
                work = {k: abs(v - ref["work_kn_m"][k]) for k, v in n["work_kn_m"].items()}
                common.append(dict(progress=n["z"][-1], factor=n["factor"], fields=row,
                                   work_difference_kn_m=work))
            parent = state
    finally:
        op.close()
    assert common
    return dict(physical_classification="D", common=common, restart_exact=True,
                paths={key:dict(accepted=len(result["points"]), last=result["last_progress"],
                               stop=result["stop"], maximum_relative_force=max(
                                   (p["relative_force"] for p in result["points"]), default=None))
                       for key, result in (("native", native), ("independent", reference))})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    p.add_argument("bridge", type=Path)
    args = p.parse_args()
    with (args.folder / "comparison.json").open("x", encoding="utf-8") as stream:
        json.dump(compare(args.folder, args.bridge), stream, indent=2, allow_nan=False)
