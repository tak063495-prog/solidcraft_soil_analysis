"""Probe fixed-parent interface equilibria during reverse and reload on three spaces."""

import inspect
import json
from pathlib import Path
import sys

import numpy as np

root = Path.cwd()
sys.path.insert(0, str(root / "interim/2026-11-25-approved-six"))
import c02_reference as ref


def probe(v):
    name = v["stage"]["name"]
    if name not in ("reverse", "reload") or v["step"] != 128:
        return
    x = v["x"]
    n = v["n"]
    evaluate = v["evaluate"]

    def jacobian(jump):
        width = 2e-9
        offsets = width * np.concatenate((np.eye(3), -np.eye(3)))
        q, _, _ = ref.project(
            v["g"].system,
            v["old_jump"],
            v["traction"],
            jump[:, None, :] + offsets,
            v["counts"],
        )
        local = ((q[:, :3] - q[:, 3:]) / (2 * width)).transpose(0, 2, 1)
        return v["condensed"] + v["c"].T @ (
            v["g"].area[:, None, None] * (local @ v["ca"])
        ).reshape(3 * n, -1)

    k = jacobian(v["jump"])
    _, s, vt = np.linalg.svd(k)
    rows = []
    for direction in (vt[-1], vt[0]):
        for sign in (-1, 1):
            y = x + sign * 1e-5 * direction
            for it in range(20):
                j, _, flow, labels, r = evaluate(y)
                if np.linalg.norm(r) / v["scale"] < 1e-10:
                    break
                dy = np.linalg.solve(jacobian(j), -r)
                for power in range(16):
                    trial = y + dy / 2**power
                    if np.linalg.norm(evaluate(trial)[-1]) < np.linalg.norm(r):
                        y = trial
                        break
                else:
                    raise ValueError("diagnostic direction did not descend")
            j, _, flow, labels, r = evaluate(y)
            rows.append(
                dict(
                    du=np.linalg.norm(x - y),
                    flow=np.max(abs(flow - v["flow"])),
                    relative=np.linalg.norm(r) / v["scale"],
                    distinct_labels=int(np.count_nonzero(labels != v["status"])),
                    enough_partitions=bool(
                        np.all(
                            np.linalg.norm(j - v["old_jump"], axis=1)
                            <= v["counts"] * v["width"]
                        )
                    ),
                )
            )
    print(
        json.dumps(
            dict(
                p=v["degree"],
                phase=name,
                from_step=v["step"] - 1,
                to_step=v["step"],
                steps=v["steps"],
                smin=s[-1],
                smax=s[0],
                sym_min=np.linalg.eigvalsh((k + k.T) / 2)[0],
                probes=rows,
            )
        ),
        flush=True,
    )


source = inspect.getsource(ref.reference)
line = "            u = fixed*factor"
assert source.count(line) == 1
ref.__dict__["probe"] = probe
exec(
    compile(
        source.replace(line, "            probe(locals())\n" + line),
        ref.__file__,
        "exec",
    ),
    ref.__dict__,
)
mesh, ground, case = [
    json.loads((root / "examples" / (name + ".json")).read_bytes())
    for name in (
        "joint-water-column-tet10-mesh",
        "ground-joint-water-mc",
        "stages-joint-water-cycle",
    )
]
for degree in (2, 3, 4):
    path = Path(sys.argv[1]) / f"roots-p{degree}"
    path.mkdir()
    cost = ref.reference(mesh, ground, case, degree, 256, 6, path)
    assert cost["accepted"] + 8 <= 2065
