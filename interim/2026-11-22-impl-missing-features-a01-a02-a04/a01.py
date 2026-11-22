"""Compare independently advanced load-controlled histories, in input GP order."""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-22-impl-t7-physical-instability"))
from parallel import Model, fixtures, refine, run, serial

def observations(mesh, case, engine, independent, model):
    layout = engine["joint_integration_points"]
    order = []
    for p in model.system.points:
        hits = [i for i, q in enumerate(layout)
                if all(p[k] == q[k] for k in ("joint_id", "face_index", "owner_a", "owner_b"))
                and np.linalg.norm(np.array(p["position_m"]) - q["position_m"]) < 1e-12
                and np.linalg.norm(np.array(p["normal"]) - q["normal"]) < 1e-12
                and abs(p["weight_m2"] - q["weight_m2"]) < 1e-14]
        if len(hits) != 1:
            raise ValueError("ambiguous oriented quadrature point")
        order.append(hits[0])
    assert len(set(order)) == len(layout)
    rows = []
    for frame in engine["frames"]:
        state = frame["state"]
        key = frame["name"], state["stage_progress"]
        if key not in independent:
            continue
        row, ref = independent[key]
        pts = [state["joints"]["points"][i] for i in order]
        mismatch = sum(p["status"] != ("open", "stick", "slip")[r]
                       for p, r in zip(pts, ref["contact_status"], strict=True))
        generations = [v["generation"] for v in state["joints"]["lifecycle"]]
        supports = all(np.array_equal(np.asarray(u)[mesh["node_sets"][c["node_set"]]] @ c["direction"],
                                      np.full(len(mesh["node_sets"][c["node_set"]]),row["factor"]*c["value_m"]))
                       for u in (state["displacement_m"],ref["u"]) for c in case["boundary"]["constraints"])
        native_increasing = any(a["stage"] == key[0] and a["accepted"]
                                and a["progress_from"] < a["progress_to"] == key[1]
                                for a in engine["attempts"])
        coordinates = [progress for stage, progress in independent if stage == key[0]]
        increasing = native_increasing and all(a < b for a, b in zip([0.]+coordinates, coordinates))
        same_factor = state["boundary_factor"] == row["factor"]
        same = increasing and same_factor and supports and np.array_equal(generations, ref["generation"]) and mismatch == 0
        fields = {"u_m": (state["displacement_m"], ref["u"]),
                  "traction_kpa": ([p["traction_kpa"] for p in pts],ref["traction"]),
                  "jump_m": ([p["jump_m"] for p in pts],ref["jump"]),
                  "stress_kpa": ([p["effective_stress_compression_positive_kpa"] for p in state["material_states"]],ref["stress"])}
        rows.append(dict(stage=key[0],progress=key[1],increasing=increasing,same_load_factor=same_factor,supports_equal=supports,
                         native_generations=generations,reference_generations=ref["generation"].tolist(),
                         activity_mismatches=mismatch,absolute_l2={k:float(np.linalg.norm(np.asarray(a)-b))
                         for k,(a,b) in fields.items()} if same else None))
    return rows

def study(folder, cli, open_fraction=None):
    folder.mkdir(parents=True,exist_ok=False)
    result = {}
    for size in (18,144):
        mesh, ground, case = fixtures.stacked(steps=8)
        if size == 144:
            mesh = refine(mesh)
        if open_fraction is not None:
            # A separate load path with its own accepted increments. The two
            # engines use identical endpoints and the original eight fractions.
            case["stages"]=case["stages"][:3]
            origin=case["stages"][1]["target_factor"]
            case["stages"][2]["target_factor"]=origin+open_fraction*(case["stages"][2]["target_factor"]-origin)
        dest = folder/str(size)
        dest.mkdir()
        inputs = [dest/(k+".json") for k in ("mesh","ground","case")]
        for p,x in zip(inputs,(mesh,ground,case),strict=True):
            p.write_text(json.dumps(x,default=serial),encoding="utf-8")
        runfile = dest/"run.json"
        with (dest/"native.log").open("w") as log:
            proc=subprocess.run([cli,"solve-stages",*inputs,runfile],stdout=log,stderr=subprocess.STDOUT,timeout=86400)
        native=json.loads(runfile.read_bytes())
        assert (proc.returncode==0 and native["complete"]) or native["failure_reason"]
        if open_fraction is None:
            assert not native["complete"]
        results={"native_stop":native["failure_reason"],"native_exit":proc.returncode,
                 "native_accepted":[[v["stage"],v["progress_to"]] for v in native["attempts"] if v["accepted"]]}
        for policy in ("centered","installed_stick"):
            model=Model(mesh,ground,case)
            states={}
            def sink(row, fields):
                states[(row["stage"],row["progress"])]=(row,fields)
            reason=None
            try:
                run(model,8,sink,birth_predictor=policy)
            except (ArithmeticError,np.linalg.LinAlgError) as e:
                reason=str(e)
            results[policy]=dict(stop=reason,accepted=list(states),
                                 comparison=observations(mesh,case,native,states,model))
        result[str(size)]=results
    return result

if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("folder",type=Path)
    p.add_argument("cli",type=Path)
    p.add_argument("--open-fraction",type=float,choices=(0.02,))
    a=p.parse_args()
    data=study(a.folder,a.cli.resolve(),a.open_fraction)
    (a.folder/"results.json").write_text(json.dumps(data,indent=2,allow_nan=False)+"\n",encoding="utf-8")
