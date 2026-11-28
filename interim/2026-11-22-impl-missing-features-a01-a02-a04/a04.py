"""Run both choices on the same wet column; advance the reference from inputs alone."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-20-impl-t1-wet-joint-richards"))
from study import cases, EX, read, write
from t1_water_balance import audit
from a04_reference import reference


def fixture():
    mesh = read(EX/"joint-water-column-mesh.json")
    ground = read(EX/"ground-joint-water-mc.json")
    ground["initialization"].update(water_level_m=None)
    ground["initialization"]["geostatic"]["hydrostatic_face_sets"] = []
    case = dict(cases())["staggered"]
    flow = case["flow"]
    day = 86400.
    for m in flow["materials"]:
        m["volume_source_per_s"] = 0.
    period = flow["analysis"]["periods"][0]
    period["end_time_s"] = 2*day
    period["boundaries"] = [
        dict(face_set=side,condition=dict(kind="head",value=dict(head_m=h)))
        for side,h in (("bottom",2.1),("top",2.0))]
    flow["settings"].update(initial_time_step_s=day, maximum_time_step_s=day)
    return mesh, ground, case, day


def compare(mesh, ground, case, run, dt):
    expected = reference(mesh,ground,case,dt,2)
    xyz = np.array(mesh["nodes"])
    faces = [tuple(sorted(tuple(p) for p in xyz[f["corners"]])) for f in run["hydraulic_faces"]]
    measurements = []
    before = next(f["state"]["fluid"] for f in reversed(run["frames"]) if f["phase"]=="initialization")
    for frame, ref in zip((f for f in run["frames"] if f["phase"]=="transient"),expected,strict=True):
        state = frame["state"]
        water = state["hydraulic"]
        assert water["time_s"] == ref["time"] and frame["staggered_iterations"] == 2
        assert all(p["status"]=="stick" for p in state["mechanical"]["joints"]["points"])
        curr = state["fluid"]
        delta = {i:v-before["joint_water_volume_m3"][i] for i,v in curr.get("joint_water_volume_m3",{}).items()}
        traces = {t["side_a"]["face"]: t for t in curr["joint_traces"]}
        assert all(t["side_b"]["face"]==i for i,t in traces.items())
        measurement=dict(time_s=ref["time"],accepted_exchanges=frame["staggered_iterations"],
                         cell_head_max_error_m=float(np.max(np.abs(water["relative_cell_head_m"]-ref["heads"]))),
                         face_head_max_error_m=float(np.max(np.abs(
                             water["relative_face_head_m"]-np.array([ref["face_heads"][k] for k in faces])))),
                         displacement_max_error_m=float(np.max(np.abs(state["mechanical"]["displacement_m"]-ref["u"]))),
                         outflow_max_error_m3_s=float(np.max(np.abs(water["cell_outflow_m3_s"]-ref["flux"]))),
                         storage_by_face=[dict(face=int(i),change_m3=v,reference_change_m3=ref["volume_change"][faces[int(i)]],
                                               both_outflows_m3_s=traces[int(i)]["side_a"]["outward_flow_m3_s"]
                                               +traces[int(i)]["side_b"]["outward_flow_m3_s"]) for i,v in delta.items()],
                         balance=frame["hydraulic_balance"])
        measurements.append(measurement)
        before=curr
    budgets={}
    for profile in ("verification","practical"):
        ledger=audit(mesh,run,profile)
        budgets[profile]=dict(passed=ledger["passed"],rows=len(ledger["rows"]),
                              interface_storage_m3=ledger["interface_storage_m3"],
                              global_rows=[r for r in ledger["rows"] if r["kind"]=="global"])
    return dict(accepted=len(measurements),end_time_s=measurements[-1]["time_s"],
                comparisons=measurements,budgets=budgets)


def study(out,cli):
    out.mkdir(parents=True,exist_ok=False)
    result={}
    runs=[]
    for name in ("zero","filled"):
        mesh,ground,case,dt=fixture()
        if name=="filled":
            case["joint_aperture_storage_m"]=0.001
        dest=out/name
        dest.mkdir()
        for key,value in (("mesh",mesh),("ground",ground),("case",case)):
            write(dest/(key+".json"),value)
        with (dest/"solve.log").open("x") as log:
            subprocess.run([cli,"solve-one-way",*[dest/(k+".json") for k in ("mesh","ground","case","run")]],
                           stdout=log,stderr=subprocess.STDOUT,check=True,timeout=86400)
        run=read(dest/"run.json")
        assert run["complete"],run["failure_reason"]
        result[name]=compare(mesh,ground,case,run,dt)
        runs.append([f["state"]["hydraulic"] for f in run["frames"] if f["phase"]=="transient"])
    result["exchange_difference"]=[dict(time_s=a["time_s"],max_cell_head_m=float(np.max(np.abs(
        np.asarray(a["relative_cell_head_m"])-b["relative_cell_head_m"]))),
        max_outflow_m3_s=float(np.max(np.abs(np.asarray(a["cell_outflow_m3_s"])-b["cell_outflow_m3_s"]))))
        for a,b in zip(*runs,strict=True)]
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("out",type=Path)
    parser.add_argument("cli",type=Path)
    args=parser.parse_args()
    write(args.out/"results.json",study(args.out,args.cli.resolve()))
