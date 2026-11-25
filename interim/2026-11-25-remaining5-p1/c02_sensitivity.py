"""Measure p, face rule and localized step-size changes at the original 448 anchors."""

import inspect
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [
    str(HERE.parent / "2026-11-25-approved-six"),
    str(HERE.parent / "2026-11-25-remaining5-diagnosis"),
]
import c02_reference as ref
from c02_series import compare, integral


# Reuse the validated equilibrium loop with an explicit time grid and a passive
# observation on common Q6 points. No observation enters equilibrium or history.
def instrument():
    src = inspect.getsource(ref.reference)
    changes = {
        "face_order,folder):": "face_order,folder,grids,observe):",
        "for step in range(1,steps+1):": "for t in grids[phase]:",
        "            t = step/steps\n": "",
        "if step%(steps//64)==0:": """observe(phase,t,g,u,status,width)
            if t*64 == int(t*64):""",
    }
    for a, b in changes.items():
        assert src.count(a) == 1
        src = src.replace(a, b)
    ns = dict(ref.__dict__)
    exec(compile(src, ref.__file__, "exec"), ns)
    return ns["reference"]


def main(out):
    out.mkdir()
    inputs = [
        json.loads((ROOT / "examples" / (name + ".json")).read_bytes())
        for name in (
            "joint-water-column-tet10-mesh",
            "ground-joint-water-mc",
            "stages-joint-water-cycle",
        )
    ]
    stages = ["geostatic"] + [s["name"] for s in inputs[2]["stages"]]
    solve = instrument()
    switch = [set() for _ in stages]
    records, costs = {}, {}
    paths = [(p, q, 1, "switch") for p, q in ((3, 6), (2, 6), (4, 6), (3, 8), (3, 10))]
    paths.extend((3, 6, sub, mode) for mode in ("switch", "reverse") for sub in (2, 4))
    for p, q, sub, mode in paths:
        tag = f"p{p}-q{q}-{mode}{sub}"
        folder = out / tag
        folder.mkdir()
        grids = [[] for _ in stages]
        for phase, grid in enumerate(grids):
            for k in range(1, 257):
                selected = (
                    k in switch[phase]
                    if mode == "switch"
                    else stages[phase] == "reverse"
                )
                n = sub if selected else 1
                grid.extend((k - 1 + j / n) / 256 for j in range(1, n + 1))
        common = ref.geometry(*inputs, p, 6)
        jump = np.zeros((len(common.area), 3))
        traction, plastic = jump.copy(), jump.copy()
        samples, labels = [], []
        last = None

        def observe(phase, t, g, u, status, width):
            nonlocal last, jump, traction, plastic
            if (p, q, sub) == (3, 6, 1) and stages[phase] in ("reverse", "reload"):
                if last is not None and np.any(status != last):
                    switch[phase].add(round(t * 256))
            last = status.copy()
            if q != 6 or (p, sub) == (3, 1):
                target = common.scalar_c @ u.reshape(-1, 3)
                counts = np.maximum(
                    1, np.ceil(np.linalg.norm(target - jump, axis=1) / width)
                ).astype(int)
                tr, flow, st = ref.project(
                    common.system, jump, traction, target[:, None, :], counts
                )
                jump, traction = target, tr[:, 0]
                plastic += flow[:, 0]
                if t * 64 == int(t * 64):
                    samples.append(plastic.copy())
                    labels.append(st[:, 0].copy())

        costs[tag] = solve(*inputs, p, 256, q, folder, grids, observe)
        with np.load(folder / "slip.npz") as raw:
            area = integral(raw)
            if (p, q, sub) == (3, 6, 1):
                np.testing.assert_allclose(samples, raw["slip"], atol=5e-15, rtol=0)
                np.testing.assert_array_equal(labels, raw["status"])
            common_fields = (
                dict(
                    slip=np.array(samples),
                    status=np.array(labels),
                    position=common.positions,
                    parent=common.parents,
                    weight=common.area,
                )
                if q != 6
                else {k: raw[k] for k in raw}
            )
        assert common_fields["slip"].shape == (448, 72, 3) and area.shape == (448, 3)
        np.savez_compressed(
            folder / "principal.npz",
            **common_fields,
            signed_area_integral=area,
            phase=np.repeat(stages, 64),
            progress=np.tile(np.arange(1, 65) / 64, 7),
        )
        records[p, q, sub, mode] = common_fields, area

    def pair(a, b):
        f, x = records[(*a, "switch") if len(a) == 3 else a]
        g, y = records[(*b, "switch") if len(b) == 3 else b]
        row = compare(f, g)
        row.pop("integral_relative")
        row["native_area_relative"] = np.linalg.norm(x - y) / np.linalg.norm(y)
        row["by_stage"] = {}
        for i, name in enumerate(stages):
            s = slice(i * 64, (i + 1) * 64)
            w = np.sqrt(g["weight"])[None, :, None]
            row["by_stage"][name] = dict(
                field_difference_norm=np.linalg.norm((f["slip"][s] - g["slip"][s]) * w),
                integral_difference_norm=np.linalg.norm(x[s] - y[s]),
            )
        return row

    report = dict(
        costs=costs,
        refined_base_intervals={name: sorted(k) for name, k in zip(stages, switch)},
        degree=[pair((p, 6, 1), (p + 1, 6, 1)) for p in (2, 3)],
        face=[pair((3, q, 1), (3, q + 2, 1)) for q in (6, 8)],
        time=[pair((3, 6, s), (3, 6, 2 * s)) for s in (1, 2)],
        reverse=[
            pair((3, 6, 1), (3, 6, 2, "reverse")),
            pair((3, 6, 2, "reverse"), (3, 6, 4, "reverse")),
        ],
    )
    (out / "summary.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]))
