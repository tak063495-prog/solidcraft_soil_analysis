"""Verify saved fold/terminal equilibria and summarize every path trial."""
import gzip
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
P2 = HERE.parent / "2026-11-25-remaining5-p2"
sys.path.insert(0,str(P2))
from nc_continue import Section, internal, prepare


def physical_check(mesh, run, cp, q, off, material):
    # Evaluate the original dimensionless ellipse after division by pc.
    # At the smallest accepted stresses, squaring in kPa first underflows.
    force = np.zeros((len(mesh["nodes"]),3))
    worst, min_p, min_pc = -np.inf, np.inf, np.inf
    laws = {m["id"]:m["mcc"] for m in run["ground"]["materials"]}
    for e in np.flatnonzero(cp["active_elements"]):
        points = material[off[e]:off[e+1]]
        stress = np.array([s["effective_stress_compression_positive_kpa"] for s in points])
        pc = np.array([s["preconsolidation_pressure_kpa"] for s in points])
        assert np.all(pc>0) and all(s["void_ratio"]>0 for s in points)
        scaled = stress/pc[:,None]
        p = scaled[:,:3].mean(axis=1)
        dev = scaled.copy()
        dev[:,:3] -= p[:,None]
        m = laws[cp["element_material_ids"][e]]["critical_state_m"]
        y = 1.5*(np.square(dev[:,:3]).sum(axis=1)+2*np.square(dev[:,3:]).sum(axis=1))/m**2 + p*(p-1)
        assert (p>0).all() and (y<=1e-8).all()
        worst = max(worst,float(y.max()))
        min_p = min(min_p,float(stress[:,:3].mean(axis=1).min()))
        min_pc = min(min_pc,float(pc.min()))
        qp = q[e]
        f = -np.einsum("qai,qa,q->i",qp["b"],stress,qp["weight"])
        force[qp["nodes"]] += f.reshape(-1,3)
    if min_pc>1e-150:
        original,_ = internal(mesh,run["ground"],material,cp["active_elements"],q,off,
                              material_ids=cp["element_material_ids"])
        np.testing.assert_allclose(force,original,rtol=0,atol=1e-12)
    return force,worst,min_p,min_pc


def audit(inputs, routes):
    m = Section(inputs,None)
    mesh = json.loads((inputs/"mesh.json").read_bytes())
    with zipfile.ZipFile(HERE.parent / "2026-11-23-nc01-coordinate-polish/evidence.zip") as archive:
        run = json.loads(archive.read("run.json"))
    cp = run["last_committed"]
    q, off = prepare(mesh,run["integration_points"])
    u0 = np.ravel(cp["displacement_m"])
    fixed = np.setdiff1d(np.arange(len(u0)),m.free)
    results = {}
    for folder in sorted(routes.iterdir()):
        if not folder.is_dir():
            continue
        meta = json.loads((folder/"summary.json").read_bytes())
        trials = [json.loads(line) for line in (folder/"trials.jsonl").open()]
        accepted = [r for r in trials if r["accepted"]]
        assert len(trials)==meta["attempts"] and len(accepted)==meta["accepted"]
        assert len(trials)+1 <= min(2000,16520)
        checks = {}
        for name in ("fold","last"):
            with gzip.open(folder/(name+".json.gz"),"rt") as f:
                state = json.load(f)
            np.testing.assert_array_equal(np.ravel(state["displacement_m"])[fixed],u0[fixed])
            material = list(cp["material_states"])
            for i,s in zip(m.point_ids,state["active_material_states"],strict=True):
                material[i] = s
            force, y, min_p, min_pc = physical_check(mesh,run,cp,q,off,material)
            load = m.load(state["load"])
            ratio = np.linalg.norm((load-force.ravel())[m.free])/max(1.,np.linalg.norm(load))
            assert ratio<=m.settings["relative_force_tolerance"]
            checks[name]=dict(load=state["load"],independent_free_force=ratio,normalized_yield=y,
                              minimum_pressure_kpa=min_p,minimum_pc_kpa=min_pc)
        lin = [max(i["linear"],i["physical_linear"]) for r in accepted for i in r["iterations"] if "linear" in i]
        assert max(lin)<=m.settings["relative_linear_tolerance"]
        minrow = min(accepted,key=lambda r:r["load"])
        results[folder.name]=dict(attempts=meta["attempts"],accepted=len(accepted),stop=meta["stop"],
            last_lambda=meta["last_lambda"],arc_m=meta["arc_m"],minimum_lambda=minrow["load"],
            minimum_at_arc_m=minrow["arc_m"],positive=sum(r["delta_lambda"]>0 for r in accepted),
            negative=sum(r["delta_lambda"]<0 for r in accepted),
            last_tangent=accepted[-1]["lambda_tangent"],
            sign_changes_of_tangent=sum(a["lambda_tangent"]*b["lambda_tangent"]<0
                                        for a,b in zip(accepted,accepted[1:])),
            lambda_at_common_arc_m={str(x):float(np.interp(x,[0.]+[r["arc_m"] for r in accepted],
                                        [meta["start_lambda"]]+[r["load"] for r in accepted]))
                                    for x in (0.1,1.,5.,7.) if x<=meta["arc_m"]},
            last_rejection=trials[-1].get("reason"),
            max_lambda=max(meta["start_lambda"],max(r["load"] for r in accepted)),
            min_second_singular=min(r["sigma_next"] for r in accepted),
            min_bordered_singular=min(r["bordered_min"] for r in accepted),
            min_abs_load_projection_kn=min(abs(r["projection_kn"]) for r in accepted),
            max_relative_force=max(r["relative_force"] for r in accepted),max_true_linear=max(lin),
            equilibrium_audit=checks)
    return results


if __name__=="__main__":
    print(json.dumps(audit(Path(sys.argv[1]),Path(sys.argv[2])),indent=2,allow_nan=False))
