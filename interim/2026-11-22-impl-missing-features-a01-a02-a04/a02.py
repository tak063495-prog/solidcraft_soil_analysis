"""Test the archived transfer and a 1D compressive construction on its tip mesh."""
import argparse
import ast
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-21-impl-t3-tip-crossing-birth"))
from study import extract, read, System


def execute(args, log):
    with log.open("xb") as stream:
        return subprocess.run(args, stdout=stream, stderr=subprocess.STDOUT, timeout=86400).returncode


def incompatible_node(message, system, mu):
    match = re.search(r"physical_dofs=(\[.*?\]) requested_kn=(\[.*?\])", message)
    rows, values = (ast.literal_eval(group) for group in match.groups())
    force = np.zeros(system.c.shape[-1])
    force[rows] = values
    # A positive side-b midside shape integrates a positive combination
    # of admissible tractions. A violating nodal force rules out every GP fit.
    for node, f in enumerate(force.reshape(-1, 3)):
        shape = system.c[:,2,3*node+2]
        if min(shape) >= 0 and max(shape) > 0 and np.linalg.norm(f[:2]) > -mu*f[2]:
            yield dict(node=node, force_kn=f.tolist(), minimum_shape=float(min(shape)),
                       shear_kn=float(np.linalg.norm(f[:2])), limit_kn=float(-mu*f[2]))


def study(out, cli, probe):
    out.mkdir(parents=True,exist_ok=False)
    extract(out)
    result = {}
    old = out/"fresh/refusals/body-across-tip-transfer_force"
    for policy in ("transfer_force", "transfer_force_face_mean"):
        results = {}
        for label in ("original", "uniaxial"):
            mesh, ground, case = [read(old/(n+".json")) for n in ("mesh","ground","case")]
            definition=case["joint_system"]["joints"][0]
            definition["birth_policy"]=policy
            if label == "uniaxial":
                # A separate manufactured column. Both the old and new GP transfer
                # use this exact input, including the existing top/bottom supports.
                for mat in ground["materials"]:
                    mat["mohr_coulomb"]["poisson"]=0.
                geo=ground["initialization"]["geostatic"]
                geo["initial_boundary_factor"]=-0.05
                geo["base_loads"]["surface_loads"][0]["pressure_kpa"]=0.
                for state in geo["seed"]["states"].values():
                    state["effective_stress_compression_positive_kpa"]=[0.]*6
                for bc in case["boundary"]["constraints"]:
                    if bc["direction"][2]==0.:
                        bc["value_m"]=0.
            dest=out/(label+"-"+policy)
            dest.mkdir()
            inputs=[dest/(n+".json") for n in ("mesh","ground","case")]
            for path,value in zip(inputs,(mesh,ground,case),strict=True):
                path.write_text(json.dumps(value),encoding="utf-8")
            log=dest/"solve.log"
            runfile=dest/"run.json"
            rc=execute([cli,"solve-stages",*inputs,runfile],log)
            system=System(mesh,case["joint_system"])
            if label=="original":
                diagnostic=log.read_text(encoding="utf-8")
                assert rc==1 and not runfile.exists() and "Coulomb cone" in diagnostic
                nodes=list(incompatible_node(diagnostic,system,definition["material"]["friction_coefficient"]))
                assert nodes
                results[label]=dict(exit_code=rc,accepted_joint_birth=0,progress=0.,
                                    infeasible_nodes=nodes)
            else:
                data=read(runfile)
                assert rc==0 and data["complete"]
                state=data["last_committed"]["joints"]
                q=np.array([p["birth"]["initial_traction_kpa"] for p in state["points"]])
                # Lower soil is already equilibrated under its own weight. After
                # placement the upper self weight adds a triangular stress diagram.
                z=np.array(mesh["nodes"])[:,2]
                upper=z.max()-system.points[0]["position_m"][2]
                top=next(bc for bc in case["boundary"]["constraints"]
                         if bc["node_set"]=="top" and bc["direction"]==[0.,0.,1.])
                material=ground["materials"][0]
                traction=(top["value_m"]*geo["initial_boundary_factor"]*material["mohr_coulomb"]["young_kpa"]
                          -material["unit_weight_kn_m3"]*upper**2/2)/np.ptp(z)
                ref=np.array([0.,0.,traction])
                force=np.einsum("qag,a,q->g",system.c,ref,system.weight)
                # Native layout and weights are read in the native order.
                actual=np.zeros_like(force).reshape(-1,3)
                for p,v in zip(data["joint_integration_points"],q,strict=True):
                    load=np.outer(p["shape"],v)*p["weight_m2"]
                    actual[p["nodes_a"]]-=load
                    actual[p["nodes_b"]]+=load
                counts={name:sum(a["accepted"] and a["stage"]==name for a in data["attempts"])
                        for name in ("body_fill","birth")}
                results[label]=dict(accepted=counts,progress=data["last_committed"]["stage_progress"],
                                    analytic_traction_kpa=ref.tolist(),gp_count=len(q),
                                    max_gp_difference_kpa=float(np.max(np.abs(q-ref))),
                                    max_nodal_force_difference_kn=float(np.max(np.abs(actual.ravel()-force))))
        patch=out/"fresh/projection"
        case=read(patch/"case.json")
        case["joint_system"]["joints"][0]["birth_policy"]=policy
        config=out/(policy+".json")
        config.write_text(json.dumps(case),encoding="utf-8")
        controls={}
        for name in ("good","tension","outside-circle","balanced-unreproducible-tip"):
            path=out/(policy+"-"+name+".json")
            rc=execute([probe,patch/"mesh.json",config,patch/(name+"-request.json"),path],path.with_suffix(".log"))
            r=read(path)
            assert r["accepted"] == (name=="good") and rc == (0 if r["accepted"] else 1)
            if not r["accepted"]:
                expected="every free nodal residual" if name=="balanced-unreproducible-tip" else "Coulomb cone"
                assert expected in r["error"]
            controls[name]=dict(accepted=r["accepted"],max_error_kn=r["result"]["maximum_force_error_kn"]
                                if r["accepted"] else None)
        # The two faces carry distinct, strictly compressive constant tractions.
        # The shared tip is one physical DOF, so its two virtual works cancel.
        model = System(read(patch/"mesh.json"),case["joint_system"])
        field = np.array([[0.,0.,-100. if p["face_index"]==0 else -10.]
                          for p in model.points])
        request = read(patch/"good-request.json")
        request["residual_kn"] = np.einsum("qag,qa,q->g",model.c,field,model.weight).tolist()
        prescribed = out/(policy+"-two-face-request.json")
        prescribed.write_text(json.dumps(request),encoding="utf-8")
        path = out/(policy+"-two-face-result.json")
        rc = execute([probe,patch/"mesh.json",config,prescribed,path],path.with_suffix(".log"))
        trial = read(path)
        enabled = policy=="transfer_force_face_mean"
        assert trial["accepted"]==enabled and rc==(0 if enabled else 1)
        if enabled:
            recovered = np.array([p["birth"]["initial_traction_kpa"] for p in trial["result"]["state"]["points"]])
            exact = [[0.,0.,-100. if p["face_index"]==0 else -10.] for p in trial["layout"]]
            controls["two_compressions"] = dict(accepted=True,
                max_gp_error_kpa=float(np.max(np.abs(recovered-exact))),
                max_error_kn=trial["result"]["maximum_force_error_kn"])
        else:
            assert "Coulomb cone" in trial["error"]
            controls["two_compressions"] = dict(accepted=False)
        results["projection_controls"]=controls
        result[policy]=results
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    for key in ("out","cli","probe"):
        parser.add_argument(key,type=Path)
    args=parser.parse_args()
    result=study(args.out,args.cli.resolve(),args.probe.resolve())
    (args.out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
