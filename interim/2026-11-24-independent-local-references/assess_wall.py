"""Apply the approved two-axis Kirsch wall criterion to the archived FEM measurements."""

import argparse
import json
from pathlib import Path


def assess(grids, criteria):
    [grid] = [
        g
        for g in grids
        if (g["theta_cells"], g["radial_cells"], g["axial_cells"]) == (24, 16, 1)
    ]
    [wall] = [s for s in grid["stages"] if s["release"] == 1]
    reference = criteria["kirsch"]["wall_hoop_compression_kpa_horizontal_vertical"]
    stresses = [point["wall_kpa"] for point in wall["extrapolation"]]
    if (
        grid["retained_tets"] != 2304
        or len(stresses) != 2
        or reference != [75000, 15000]
        or wall["wall_reference_kpa"] != reference
    ):
        raise ValueError("unexpected Kirsch case")
    errors = [
        (actual - ref) / ref for actual, ref in zip(stresses, reference, strict=True)
    ]
    return dict(
        reference_kpa=reference,
        wall_kpa=stresses,
        signed_relative_errors=errors,
        relative_limit_inclusive=0.03,
        passed=all(abs(e) <= 0.03 for e in errors),
        decision="ユーザー決定 2026-11-24 14:10–14:11 JST。問題があれば後から変更できる。",
        analytical_source=criteria["kirsch"]["publication"],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    interim = Path(__file__).resolve().parent.parent
    data = json.loads(
        (interim / "2026-11-24-approved-four/kirsch-wall.json").read_bytes()
    )
    spec = json.loads(
        (interim / "2026-11-24-t9d-benchmarks/criteria.json").read_bytes()
    )
    args.output.write_text(
        json.dumps(assess(data, spec), indent=2, ensure_ascii=False) + "\n",
        encoding="utf8",
    )
