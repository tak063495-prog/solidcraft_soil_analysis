"""Read-only extraction from hash-pinned legacy JSON/text. No solver imports.

Original h-series stress means FULL tensor norm, not one depth component.
Original 'oscillating monitor' is static construction uz across spatial meshes.
"""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = "3b1133d3be9905dc24f62c27b5b1eb456a83f2e2"


def collect(root=ROOT):
    manifest = json.loads((HERE/"legacy-sources.json").read_text(encoding="utf-8"))
    verified = {}
    for source in manifest["sources"]:
        raw = (root/source["path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != source["sha256"]:
            raise ValueError("legacy source missing/changed: "+source["path"])
        verified[source["path"]] = raw.decode("utf-8")
    def old(name):
        path = next(p for p in verified if p.endswith("/"+name))
        return json.loads(verified[path]),path
    diagnostic,diag_path = old("plane-diagnostic.json")
    plane,plane_path = old("plane.json")
    finite,finite_path = old("depth-finite.json")
    fine,fine_path = old("depth-finite-xz4.json")
    def field(report):
        return next(r for r in report["field_difference_orders"] if r["stage"]=="cut" and r["field"]=="stress")
    ps,fs,fs4 = field(plane),field(finite),field(fine)
    mon = diagnostic["monitor"]
    common = {"geometry":"1 m cube; x/z in-plane, y is extrusion depth, z up",
              "material":{"E_kPa":10000.,"nu":.25,"gamma_kN_m3":20.,"c_kPa":2.,
                          "phi_deg":20.,"psi_deg":0.,"flow":"non_associated_MC","Rankine":"none"},
              "initial":"zero effective stress and plastic history, void ratio .7; self-weight equilibration",
              "boundary":"x=0/1 ux=0, z=0 uz=0; y constraint differs per row",
              "load":"self-weight then Cut body removal, remaining face pressure 0 kPa; no surcharge/water",
              "time":"static construction progress; physical seconds/mass/damping 未定義",
              "increment":1/16,"minimum_increment":1/4096,"material_substep":1e-4,
              "half_control":{"increment":1/32,"material_substep":5e-5},
              "source":"tools/plastic_spatial_convergence.py:74-125; docs/PLASTIC_SPATIAL_CONVERGENCE.ja.md"}
    gates = {"relative_force_residual_max":1e-8,"global_force_error_max":1e-8,
             "relative_linear_tolerance":1e-10,"normalized_yield_max":1e-7,"engineering_strain_max":.05,
             "field_extrapolation":"Dfine < 0.7*Dcoarse, both positive, p=log2(Dcoarse/Dfine) in [0.5,8]",
             "scalar_extrapolation":"same-sign nonzero differences, |b|<.7|a|, p in [.5,8]; variation above 1e-8*max|values|",
             "source":"tools/plastic_spatial_convergence.py:38-65,272-274,294-301",
             "absolute_stress_accuracy_tolerance":None,
             "new_elastic_dynamic_reference_tolerance":None}
    rows = [
        {"id":"extrusion-stress","state":"判定不能 / 未達","conditions":{
            "depth_boundary":"uy=0 everywhere","cut":"x<.5,z>.5 all y; active volume .75 m3",
            "mesh":"in-plane x/z=2/4/8; y=1 fixed; 24/96/384 Tet4"},
         "comparison":"cut full 6-component tensor Frobenius RMS over common active volume; shear squared weight 2",
         "values_kPa":ps["difference_norms"],"difference_ratio":ps["difference_norms"][1]/ps["difference_norms"][0],
         "original_order":ps["order"],"original_estimate":ps["finest_error_norm_estimate"],
         "component_diagnostic":diagnostic["stress"],
         "reason":"ratio >= .7; an individual component improvement cannot replace full tensor gate",
         "source_pointers":[plane_path+"#/field_difference_orders/5",diag_path+"#/stress"],
         "same_condition_independent_MC_reference":"未確認; elastic Lamé benchmark is a different problem",
         "accepted":None},
        {"id":"oscillating-static-monitor","state":"判定不能 / 未達","conditions":{
            "same_as":"extrusion-stress","position_m":mon["position_m"],
            "physical_time_s":"未定義","inertial_mass":"未定義","damping":"未定義"},
         "comparison":"static cut uz minus initial-equilibrium uz at same material point, across x/z meshes",
         "values_m":mon["construction_uz_m"],"signed_changes_m":mon["construction_signed_changes_m"],
         "initial_uz_m":mon["initial_uz_m"],"total_uz_m":mon["total_uz_m"],
         "half_control_values_m":mon["half_control_richardson"]["values"],
         "original_order":mon["original_richardson"]["order"],"original_estimate":mon["original_richardson"]["estimate"],
         "reason":mon["original_richardson"]["reason"],
         "velocity_m_s":"未定義: no physical time","acceleration_m_s2":"未定義: no physical time",
         "period_s":"未定義: nonmonotone spatial sequence is not a dynamic oscillation",
         "source_pointers":[diag_path+"#/monitor"],
         "same_condition_independent_dynamic_reference":"未定義; SDOF/rod are separate explicit dynamic benchmarks",
         "accepted":None},
        {"id":"finite-width-depth-stress","state":"判定不能 / 未達","conditions":{
            "depth_boundary":"only y=0/1 uy=0, interior uy free; NOT traction-free ends",
            "cut":"x<.5,z>.5,y<.5 fixed; active volume .875 m3",
            "series_1":"x/z=2 fixed, y=2/4/8; 48/96/192 Tet4",
            "series_2":"x/z=4 fixed, y=2/4/8; 192/384/768 Tet4; do not mix"},
         "comparison":"full tensor Frobenius RMS after depth MESH refinement at FIXED physical 1 m depth; not depth-component-only or physical-length sensitivity",
         "series_1_kPa":fs["difference_norms"],"series_1_ratio":fs["difference_norms"][1]/fs["difference_norms"][0],
         "series_2_kPa":fs4["difference_norms"],"series_2_ratio":fs4["difference_norms"][1]/fs4["difference_norms"][0],
         "original_orders":[fs["order"],fs4["order"]],
         "reason":"both ratios >= .7; fixed coarse in-plane mesh also prevents general 3D convergence claim",
         "source_pointers":[finite_path+"#/field_difference_orders/5",fine_path+"#/field_difference_orders/5"],
         "same_condition_independent_MC_reference":"未確認; finite elastic box has different material/load/end constraints",
         "accepted":None}]
    # New diagnostic of immutable numbers, not a new solve or extrapolation.
    if any(r["accepted"] is not None for r in rows) or ps["order"] is not None or fs["order"] is not None or fs4["order"] is not None:
        raise ValueError("legacy indeterminate state unexpectedly changed")
    if mon["construction_signed_changes_m"][0]*mon["construction_signed_changes_m"][1] >= 0:
        raise ValueError("legacy monitor sign reversal changed")
    return {"schema":"p6-03-legacy-rows/1","base_commit":BASE,"p6_03_complete":False,"accepted":None,
            "sources_verified":manifest["sources"],"common_conditions":common,"original_gates":gates,"rows":rows,
            "axis_warning":"Legacy y is depth; new reference z is depth. Legacy stresses compression-positive kPa; new stresses tension-positive Pa. No automatic conversion.",
            "ledger_30":"ledger-30-extract.json; TODO:P6-03 differs from STATUS:P6-03"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output",type=Path)
    args = p.parse_args()
    result = collect()
    with args.output.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
        f.write("\n")


if __name__ == "__main__":
    main()
