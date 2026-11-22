"""Scoped T8 decisions from the small studies. Original rows always retain their status."""
import argparse
import json
from pathlib import Path

from refinement import richardson
from branches import evaluate
from sealed import HERE

T4 = HERE.parent / "2026-11-21-impl-t4-offset-plate-composite"
T6 = HERE.parent / "2026-11-21-impl-t6-offset-kink-joint"
THIN = HERE.parent / "2027-01-14-b02-fixed-offset-thickness-declaration-order"


def read(path):
    return json.loads(path.read_bytes())


def limited(scope, check):
    return dict(scope=scope, status="限定受入" if check else "未達")


def withheld(scope):
    return dict(scope=scope, status="判定保留（ユーザー判断待ち）")


def item(option, *controls, cap=None):
    return dict(option=option, original="未達", limited=list(controls), scale_limit=cap)


def assess(work):
    a1=evaluate()
    a2=read(work/"a02/results.json")
    a3=read(work/"a03/results.json")
    a4=read(work/"a04/budget-summary.json")
    water=read(work/"t8-water/affine.json")
    patch=read(work/"b02-patch/patch-results.json")
    flat=read(work/"b01/attached/measurements.json")
    modes=read(work/"b03/controls.json")
    shell=read(work/"b01/standard/measurements.json")
    thin=read(THIN/"measurements.json")
    c1=read(work/"c01-results.json")
    c2=read(work/"c02-extended.json")
    c3=read(work/"c03/results.json")
    restart=read(work/"c03/matrix/restarts/verification.json")
    c4=read(work/"c04-results.json")
    nc=read(next(HERE.parent.glob("*t9*/validation.json")))
    if nc["physical_classification"] != "D" or nc["nc01_complete"]:
        raise ValueError("NC-01 source classification changed")
    a2_transfer = (a2["projection/good"]["accepted"] and a2["projection/good"]["max_all_nodal_error_kn"] <= 1e-10
                   and all(not a2["projection/"+k]["accepted"] for k in ("outside-circle","tension","balanced-unreproducible-tip"))
                   and all(v["virtual_work_difference_kn_m"] <= 1e-10 for v in a2["projection_controls"].values())
                   and a2["whole_faces"]["transfer_force"]["existing_gates_passed"])
    placed = read(work/"t8-water/placement.json")
    # The closed-form tet/face unit controls have their original roundoff checks
    # in tests/t2_soil_birth.rs. Expose the moved-plane observations separately.
    placement_errors = {key:max(abs(row["difference"][key]) for row in placed)
                        for key in placed[0]["difference"]}
    water_ok = len(water)==6 and all(
        all(p["passed"] is True for p in row["water"]["profiles"].values()) and
        all(abs(r["pressure_error_kpa"])<=1e-12 and r["max_nodal_force_error_kn"]<=1e-10 and
            max(abs(r["side_a_head_error_m"]), abs(r["side_b_head_error_m"]))<=1e-10 and
            max(abs(r["a_flow_error_m3_s"]),abs(r["b_flow_error_m3_s"]))<=1e-14 for r in row["traces"])
        for row in water)
    dry_frames = [f for f in patch["frames"] if "local_nodal_work_difference" in f]
    legacy=read(T4/"measurements.json")["b02"]
    original_joins=read(T6/"b03-replay.json")
    items = {
        "A-01":item("②", dict(scope="18/144 Tet4、BIRTH/open、同一進捗・支持・世代・活動枝のみ",status="限定確認"),
                    withheld("平行2面の停止荷重・全場・区間仕事")),
        "A-02":item("①＋②限定",limited("24 Tet10、適合小面Aの許容力・全節点力・仮想仕事",a2_transfer),
                    limited("S1-B 無応力新設・S2 部分再生、独立全場9/17状態",
                            all(a2[k]["limited_field_gates_passed"] for k in ("body-across-tip-stress_free","partial-face-rebirth")))),
        "A-03":item("②",limited("12 Tet10 定常明示配置・12 Tet4 搬入/搬出・異世代再開",
                               all(a3[k]["complete"] for k in ("two_region","post_darcy","saturated_mcc","lifecycle"))
                               and all(v["passed"] for v in a3["lifecycle"]["budget"]["profiles"].values())
                               and all(v["identical"] for v in a3["lifecycle"]["checkpoint_restarts"])),
                    withheld("定常処方の全域積算水量")),
        "A-04":item("①＋②",limited("12 Tet4・終点0.25 s、同一水頭／非ゼロ横断流・3刻みの各面荷重と流量",water_ok),
                    limited("T1小型6経路の水収支、ゼロ逆伝達2反復",
                            len(a4)==6 and all(v["passed"] for x in a4.values() for v in x["profiles"].values())
                            and set(a4["staggered"]["iterations"])=={2})),
        "B-01":item("②",limited("12 Tet10・2 Tri6 の乾燥付着・凍結・別プロセス再開",
                               flat["tie_max_m"]<=1e-10 and flat["frozen_equal"] and all(flat["byte_equal"])),
                    withheld("8/32/128枚の曲面変位・断面量・最低次数")),
        "B-02":item("①＋②",limited("12 Tet10・T4乾燥パッチの既存代数条件と再開",
                               len(dry_frames)>0 and max(abs(f["local_nodal_work_difference"]) for f in dry_frames)<1e-14
                               and patch["altered_archive_rejected"] and all(r["identical"] for r in patch["restarts"])),
                    dict(scope="12/96の共通物理点と固定面積のN/M/Q・仕事",status="限定確認"),
                    cap="旧768の全場・再組立を新実行版で確認不可・未達"),
        "B-03":item("①＋②限定",limited("旧共有節点の腕整合・既存21経路",original_joins["limited_gates"]==21),
                    limited("構造単独、支持固定・剛結/ヒンジ/正ばねの27対照の零モード",
                            len(modes)==27 and all(r["modes"]["free"]["nullity"]=={"rigid":6,"hinge":7,"spring":6}[r["law"]]
                            and r["modes"]["supported"]["nullity"]==0 for r in modes)),
                    withheld("実断面・同じ境界条件の公刊参照")),
        "C-01":item("②",dict(scope="12 Tet10・32GP・7工程、切替二分4・N=128/256/512、両実装",status="限定確認"),
                    withheld("全履歴の収束とGCI")),
        "C-02":item("②",dict(scope="12/96 Tet10・7工程、共通160点・固定32/64/128、3写像、面積積分",status="限定確認"),
                    cap="元の3格子全場は768/6144を必要とし確認不可・未達"),
        "C-03":item("②",limited("12 Tet4・MC全7工程32増分、7地点の別PID再開",
                               restart["complete"] and restart["real_files_separate_process"] and len(restart["checks"])==7),
                    withheld("材料・接触・全体刻みを変えた全内部変数の収束")),
        "C-04":item("②",dict(scope="12 Tet4・6GP、規定水位と力学刻み別の全13工程",status="限定確認"),
                    withheld("厳密カテゴリ・履歴仕事・散逸の受入幅")),
        "NC-01":item("①",dict(scope="5 kPa全撤去・分類D・T9の判断待ち",status="未達"))
    }
    opts = dict(similar=True, complete=True, reference_converged=True,
                other_axes_fixed=True, smooth=True)
    euler = read(HERE.parent/"2026-11-22-impl-t7-physical-instability/results.json")["references"]["euler"]
    gci = {
        "Euler analytic control":richardson([r["Pcr"] for r in euler],
                                           [1/r["elements"] for r in euler], **opts),
        "B-03 regular compliance":richardson(
            [r["compliance"] for r in read(work/"b03/series.json") if not r["distorted"]],
            [1/2,1/4,1/8], **(opts | {"reference_converged":False})),
        "C-02 two meshes":dict(holds=["two grids", "traction direction reversal in old three-level record"],
                               observed_order=None, gci_fine=None, extrapolated=None),
    }
    signed = c1["directions"]
    for engine, fields in signed.items():
        for name, r in fields.items():
            holds = []
            if r["cosine"] is None or r["cosine"] <= 0:
                holds.append("direction reversal or undefined direction")
            if r["observed_order"] is None or r["observed_order"] <= 0:
                holds.append("non-positive or undefined order")
            gci[f"C-01 {engine} {name}"] = dict(holds=holds+["unresolved finite reference"],
                                                  extrapolated=None, gci_fine=None)
    gci["C-03 fixed material"] = dict(holds=["contact refinement differences increase again"],
                                       extrapolated=None, gci_fine=None)
    old_768 = read(HERE.parent/"2026-12-09-b02-tet768-reassembly-residual/measurements.json")["historical_fixed"]
    sensitivity = {
        "B-01 hypothetical 3 percent":[dict(case=r["case"], distorted=r["distorted"],
             within_example=None if r.get("relative_difference") is None else abs(r["relative_difference"])<=.03)
             for r in shell if r["n"]==8],
        "B-02 hypothetical 4e-10": {case:{s:values["max_relative"]<=4e-10
             for s,values in row["residuals"].items()} for case,row in thin["cases"].items()},
        "B-02 saved 768 at 4e-10":{k:v<=4e-10 for k,v in old_768.items() if k.endswith("reassembly_relative")},
        "C-04 hypothetical 3 percent":[dict(mechanical_steps=r["mechanical_steps"],water_steps=r["water_steps"],
             within_example=r["quantities"]["plastic_dissipation"]["history"]["relative_l2"]<=.03)
             for r in c4["measured"]]
    }
    return dict(base="120a66520b628e1a756e719c374f28d16d9e7581", passed=False, p6_complete=False,
                items=items), dict(
        branches=a1, placement_maximum_errors=placement_errors, gci=gci, sensitivity_only=sensitivity,
        roof_cylinder_hemisphere=[dict(case=r["case"],distorted=r["distorted"], signed_u=r.get("displacement"),
                                       signed_relative=r.get("relative_difference"),diagnostic=r["diagnostic"])
                                  for r in shell if r["n"]==8],
        thin_original={k:v["residuals"] for k,v in thin["cases"].items()},
        b02_saved_gate_count=legacy["bounded_field_paths_within_original_gates"],
        c01_directions=c1["directions"], c02_common_points=c2["common_points"],
        c03_contact_sensitivity={s:[c3["sensitivity"][f"m0c{i}-{s}--m0c{i+1}-{s}"]["all_fields"]["plastic_slip"]
                                    for i in range(3)] for s in ("native","independent")},
        c04_preserved_counts=c4["preserved_mismatch_counts"],
        c04_new_counts=[r["strict_mismatches"] for r in c4["measured"]],nc01=nc)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    decision, metrics=assess(args.work)
    for name, value in (("acceptance",decision),("comparisons",metrics)):
        with (args.output/(name+".json")).open("x",encoding="utf-8",newline="\n") as f:
            json.dump(value,f,indent=2,ensure_ascii=False,allow_nan=False)
            f.write("\n")
