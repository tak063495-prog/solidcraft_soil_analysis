"""Separate measured finite benchmarks from quantitative acceptance."""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CRITERIA = HERE.parent / "2026-11-24-t9d-benchmarks" / "criteria.json"


def assess():
    spec = json.loads(CRITERIA.read_bytes())
    rows = []

    def row(key, status, scope, evidence):
        rows.append(dict(id=key, status=status, scope=scope, evidence=evidence))

    def measured(key, ok, scope, evidence):
        row(key, "達成" if ok else "未達", scope, evidence)

    fields = json.loads((HERE / "c01-fields.json").read_bytes())
    for side in ("native", "independent"):
        source = f"c01-fields.json#/refinements/{side}-friction-events"
        case = fields["refinements"][f"{side}-friction-events"]
        complete = case["complete"] and case["common_times"] == 7 * 256
        measured(f"C-01:{side}:measurement", complete, "12 Tet10、3水準、7工程の共通座標での全場比較", source)
        ok = []
        for f in case["fields"].values():
            m = f["common"]
            ok.append(m["trend"] in ("identical", "decreasing_aligned") if "trend" in m
                      else m["coarse_mismatches"] == m["fine_mismatches"] == 0)
        measured(f"C-01:{side}:finite-trends", complete and all(ok),
                 "この3水準の連続20場の差分次数・方向と分類2場の一致", source)
    row("C-01:overall", "未達", "漸近域を含む全場の刻み収束", "c01-fields.json")

    band = spec["arch"]["first_peak_load_reading"]
    for n, case in json.loads((HERE / "arch.json").read_bytes()).items():
        beam = case["beam"]
        source = f"arch.json#/{n}/beam"
        fold = beam["reversal"] and beam["fold"] == "simple fold"
        measured(f"T7:{n}:reversal", beam["reversal"], "隣接する受理荷重3点での反転", source)
        measured(f"T7:{n}:fold", fold, "単純折返しの同定", source)
        measured(f"T7:{n}:figure7", fold and band[0] <= beam["fold_load"] <= band[1],
                 "単純折返し荷重が公開図の区間 885〜905 に入る", source)
    row("T7:full-curve", "未達", "全曲線の定量受入", "arch.json")

    grids = json.loads((HERE / "kirsch-wall.json").read_bytes())
    checks = [s for grid in grids for s in [grid["initial"], *grid["stages"]]]
    measured("Kirsch:equilibrium", len(checks) == 16 and all(
        0 <= s[k] <= spec["linear_relative_max"] for s in checks for k in ("free_relative", "linear_relative")
    ), "4格子の初期体と3解放段階の自由力・線形残差", "kirsch-wall.json")
    row("Kirsch:field", "未達", "孔壁・近傍場の精度についての定量受入", "kirsch-wall.json")

    paths = json.loads((HERE / "mcc-constant-g.json").read_bytes())
    expected_paths = {(n, h, ocr) for n in (12, 96) for h in (1e-4, 5e-5) for ocr in (1.6, 8)}
    measured("MCC:constant-G:endpoints", len(paths) == len(expected_paths)
             and {(r["elements"], r["material_substep"], r["ocr"]) for r in paths} == expected_paths
             and all(r["constant_shear_kpa"] == 250 and len(r["history"]) == 512
                     and all(abs(v / ref - 1) < spec["mcc"]["endpoint_relative_max_exclusive"]
                             for v, ref in zip(r["endpoint_p_q_v"], r["critical_p_q_v"], strict=True))
                     for r in paths), "8ケースの終端 p′・q・v のみ", "mcc-constant-g.json")
    row("MCC:constant-G:full-history", spec["missing_tolerance_state"],
        "軸ひずみに対する全履歴の定量受入", "mcc-constant-g.json")
    local = json.loads((HERE / "c03-local.json").read_bytes())
    row(local["id"], local["status"], local["scope"], "c03-local.json")
    row("C-03:coupled-accuracy", spec["missing_tolerance_state"], "湿潤経路の全23場と塑性仕事", local["coupled_evidence"])
    return rows


if __name__ == "__main__":
    print(json.dumps(assess(), ensure_ascii=False, indent=2))
