"""Produce the T9(d) scope ledger from numeric evidence and preregistered criteria."""
import argparse
import hashlib
import json
from pathlib import Path


def assess(folder):
    def read(name):
        return json.loads((folder/name).read_bytes())

    criteria = read("criteria.json")
    arch = read("arch-results.json")
    kirsch = read("kirsch-results.json")
    mcc = read("mcc-results.json")
    limit = criteria["linear_relative_max"]
    missing = criteria["missing_tolerance_state"]
    rows = []

    def row(name, passed, source, reason):
        rows.append(dict(id=name, state="達成" if passed else "未達", passed=passed,
                         evidence=source, reason=reason))

    def unavailable(name, source, reason):
        row(name, False, source, missing+"。"+reason)

    low, high = criteria["arch"]["first_peak_load_reading"]
    for n in (12,60,96):
        for method in ("beam","shooting"):
            case = arch[str(n)][method]
            loads = [p["load"] for p in case["curve"]]
            i = loads.index(max(loads))
            reversal = 0 < i < len(loads)-1 and loads[i-1] < loads[i] > loads[i+1]
            src = f"arch-results.json#/{n}/{method}"
            fold = reversal and case["fold"]=="simple fold"
            row(f"T9d:arch:{n}:{method}:fold", fold, src, "受理平衡列の増減と T7 の単純 fold 診断")
            peak = case["fold_load"]
            row(f"T9d:arch:{n}:{method}:figure7", fold and peak is not None and low<=peak<=high,
                src,"最初の単純 fold 荷重と公刊 Fig. 7 の読取り区間の照合")
    unavailable("T9d:arch:full-curve","criteria.json#/arch/full_curve_tolerance",
                "全曲線・頂点変位の点別許容幅は未確認。")
    if sorted(r["retained_tets"] for r in kirsch)!=[12,96]:
        raise ValueError("missing Kirsch mesh")
    for r in kirsch:
        n = r["retained_tets"]
        stages = r["stages"]
        if [s["release"] for s in stages]!=criteria["kirsch"]["release"]:
            raise ValueError("incomplete excavation history")
        checks = [r["initial"],*stages]
        ok = all(0<=s[field]<=limit for s in checks for field in ("linear_relative","free_relative"))
        row(f"T9d:kirsch:{n}:stage-equilibrium",ok,"kirsch-results.json",
            "無孔、実応力の引継ぎ、半解放、全解放の各自由力残差と線形残差")
    unavailable("T9d:kirsch:field","criteria.json#/kirsch/field_relative_tolerance",
                "解析応力・変位は既知。孔壁回復応力・変位場のFE許容幅は未確認。")
    expected = {(n,ocr,h) for n in (12,96) for ocr in (1.6,8.) for h in (1e-4,5e-5)}
    if len(mcc)!=len(expected) or {(r["elements"],r["ocr"],r["material_substep"]) for r in mcc}!=expected:
        raise ValueError("missing or repeated triaxial path")
    for ocr in (1.6,8.):
        cases = [r for r in mcc if r["ocr"]==ocr]
        errors = [abs(value/ref-1) for r in cases for value,ref in
                  zip(r["endpoint_p_q_v"],r["critical_p_q_v"],strict=True)]
        row(f"T9d:mcc:OCR{ocr:g}:critical",
            len(errors)==12 and all(e<criteria["mcc"]["endpoint_relative_max_exclusive"] for e in errors),
            "mcc-results.json","両格子・両材料内刻みで p′, q, v の各終端を解析解と比較")
        row(f"T9d:mcc:OCR{ocr:g}:free-equilibrium",
            all(0<=r["max_free_relative"]<=limit for r in cases),
            "mcc-results.json","圧密・除荷・排水軸圧縮の全増分、自由節点で再組立した力の釣合い")
    unavailable("T9d:mcc:history-accuracy","criteria.json#/mcc/history_pc_or_volume_tolerance",
                "圧密・除荷の p′c, v を持ち越し誤差を測定。公刊の前史の数値許容幅は未確認。")
    unavailable("T9d:mcc:strain-curve","criteria.json#/mcc/full_stress_strain_tolerance",
                "固定νと一定Gでは同じひずみ曲線が参照にならない。")
    row("NC-01:overall",False,"../2026-11-23-nc01-t9d-verification-proposal/ledger.json",
        "元の5 kPa 3D MCC 全撤去は未達・分類D。")
    files = ["criteria.json","arch-results.json","kirsch-results.json","mcc-results.json"]
    return {"schema":"t9d-scope-ledger/1","decision_jst":criteria["decision_jst"],
            "sha256":{name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in files},
            "rows":rows}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder",type=Path)
    args = parser.parse_args()
    print(json.dumps(assess(args.folder),indent=2,ensure_ascii=False,allow_nan=False))
