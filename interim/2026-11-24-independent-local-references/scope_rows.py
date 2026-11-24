"""Produce scoped ledger rows directly from the three verification reports."""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rows():
    wall, mcc, mc = (
        json.loads((HERE / f).read_bytes())
        for f in ("wall.json", "mcc.json", "c03.json")
    )
    output = []

    def record(key, ok, scope, evidence):
        output.append(
            dict(
                id=key, status="達成" if ok else "未達", scope=scope, evidence=evidence
            )
        )

    record(
        "Kirsch:wall:complete-release:24x16",
        wall["passed"],
        "水平・鉛直軸の壁面外挿値。解析解75／15 MPa、各3%以内。" + wall["decision"],
        "wall.json",
    )
    record(
        "Kirsch:field",
        False,
        "近傍応力・変位場の定量受入",
        "../2026-11-24-approved-four/scope-rows.json",
    )
    for i, case in enumerate(mcc["cases"]):
        name = f"MCC:{case['pc_kpa']:g}:{case['elements']}:{case['material_substep']:g}"
        source = f"mcc.json#/cases/{i}"
        for field, difference in case["maxima"].items():
            ref = mcc["refinement"][str(int(case["pc_kpa"]))][-1][field]["relative"]
            record(
                name + ":" + field,
                difference["relative"] < mcc["limit_exclusive"]
                and ref < mcc["limit_exclusive"] / 100,
                "独立実装による検証、512点の " + field,
                source,
            )
        record(
            name + ":history", case["passed"], "独立実装による検証、6量の全履歴", source
        )
    record(
        "MCC:published:history",
        False,
        "論文数値との全履歴比較。参照なし（全文待ち）",
        "../2026-11-24-approved-four/scope-rows.json",
    )
    for route, points in mc["routes"].items():
        for i, p in enumerate(points):
            key = f"C-03:{route}:GP{p['gp']}:"
            source = f"c03.json#/routes/{route}/{i}"
            record(
                key + "classification",
                p["classification_passed"],
                "独立実装による検証、全224受理点の分類",
                source,
            )
            for label, passed in zip(
                ("direct", "volumetric", "shear"), p["work_passed"], strict=True
            ):
                record(
                    key + label,
                    passed,
                    "独立実装による検証、単一材料点の累積仕事 " + label,
                    source,
                )
    record(
        "C-03:coupled-and-contact",
        False,
        "湿潤経路の全23場・全体仕事・接触遷移の定量受入。参照なし（全文待ち）",
        "../2026-11-24-approved-four/c03-local.json",
    )
    return output


if __name__ == "__main__":
    print(json.dumps(rows(), ensure_ascii=False, indent=2, allow_nan=False))
