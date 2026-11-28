"""Compare the zero-increment NC checkpoint operator without advancing the load."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from verify_high_order import internal, physical, prepare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("archive", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    args.work.mkdir()
    start = time.monotonic()
    measured = json.loads((ROOT/"interim/2026-11-23-nc01-coordinate-polish/measurements.json").read_bytes())
    for name, sha in measured["input_sha256"].items():
        assert hashlib.sha256((args.inputs/name).read_bytes()).hexdigest() == sha
    with zipfile.ZipFile(args.archive) as archive:
        raw = archive.read("run.json")
    assert hashlib.sha256(raw).hexdigest() == measured["run_sha256"]
    run_path = args.work/"archived-run.json"
    run_path.write_bytes(raw)
    run = json.loads(raw)
    mesh = json.loads((args.inputs/"mesh.json").read_bytes())
    ground, case, state = run["ground"], run["case"], run["last_committed"]
    lam = state["stage_progress"]
    assert lam == 0.5647809505462646 and not run["complete"]
    layout = run["integration_points"]
    quadrature, offsets = prepare(mesh, layout)
    live = np.array(state["active_elements"])
    live_nodes = np.unique(np.array(mesh["tetrahedra10"])[live])
    fixed = np.array(mesh["node_sets"]["bottom"])
    free = np.concatenate([np.arange(3*i,3*i+3) for i in np.setdiff1d(live_nodes, fixed)])
    assert len(free) == 348
    assert not any(live[p["owner_a"]] and live[p["owner_b"]] for p in run["joint_integration_points"])
    directions = np.eye(len(free))
    direction_file = args.work/"directions.json"
    direction_file.write_text(json.dumps(directions.tolist()))
    response_path = args.work/"native.json"
    subprocess.run([ROOT/"target/release/examples/nc_checkpoint_response",args.inputs/"mesh.json",
                    run_path,direction_file,response_path],check=True)
    native = json.loads(response_path.read_bytes())
    # The checkpoint lists the production equation ordering. Independently
    # constrained physical DOFs must match that ordering up to permutation.
    order = [(g,row[0][0]) for g,row in enumerate(native["dof_rows"]) if row]
    assert all(len(row)==1 and row[0][1]==1 for row in native["dof_rows"] if row)
    assert sorted(g for g,i in order)==list(free) and sorted(i for g,i in order)==list(range(len(free)))
    global_ids = [g for g,i in sorted(order,key=lambda pair:pair[1])]
    force, _ = internal(mesh,ground,state["material_states"],live,quadrature,offsets,
                        material_ids=state["element_material_ids"])
    final_load = physical(mesh,ground,live,state["surface_pressure_kpa"],case["pressure_on_original_surface_only"])[0].ravel()
    holding = np.asarray(run["excavation_audits"][0]["initial_holding_load_kn"])
    external = final_load+(1-lam)*holding
    residual = (external-force.ravel())[global_ids]
    nd = len(external)
    k = np.zeros((nd,nd))
    elastic_k = np.zeros_like(k)
    stencils_path = args.work/"local-stresses.json"
    subprocess.run([ROOT/"target/release/examples/nc_checkpoint_response","--local-stress-stencil",
                    run_path,stencils_path],check=True)
    stencils = {row["point"]:row for row in json.loads(stencils_path.read_bytes())}
    h = case["settings"]["tangent_perturbation"]
    materials = {m["id"]:m["mcc"] for m in ground["materials"]}
    for i in np.flatnonzero(live):
        q = quadrature[i]
        mat = materials[state["element_material_ids"][i]]
        states = state["material_states"][offsets[i]:offsets[i+1]]
        p = np.array([np.mean(s["effective_stress_compression_positive_kpa"][:3]) for s in states])
        bulk = (1+np.array([s["void_ratio"] for s in states]))*p/mat["kappa"]
        shear = 1.5*bulk*(1-2*mat["poisson"])/(1+mat["poisson"])
        # The zero-increment material map takes its elastic predictor branch.
        tensor = np.zeros((len(p),6,6))
        tensor[:,:3,:3] = (bulk-2*shear/3)[:,None,None]
        for axis in range(3):
            tensor[:,axis,axis] += 2*shear
            tensor[:,axis+3,axis+3] = shear
        idx = (3*q["nodes"][:,None]+np.arange(3)).ravel()
        elastic_k[np.ix_(idx,idx)] += np.einsum("qai,qab,qbj,q->ij",q["b"],tensor,q["b"],q["weight"],optimize=True)
        stress = np.array([stencils[j]["stress"] for j in range(offsets[i],offsets[i+1])])
        tensor = ((stress[:,1::2]-stress[:,::2])/(2*h)).transpose(0,2,1)
        k[np.ix_(idx,idx)] += np.einsum("qai,qab,qbj,q->ij",q["b"],tensor,q["b"],q["weight"],optimize=True)
    matrix = k[np.ix_(global_ids,global_ids)]
    native_k = np.asarray(native["actions_kn"]).T
    native_r = np.asarray(native["reduced_residual_kn"])
    relative_delta = np.linalg.norm(matrix-native_k)/np.linalg.norm(native_k)
    singular = np.linalg.svd(matrix,compute_uv=False)
    report = dict(checkpoint_lambda=lam,active_elements=int(live.sum()),material_points=len(layout),
                  joint_points=len(run["joint_integration_points"]),active_joint_points=0,
                  dofs=len(free),tangent_columns=len(directions),equilibrium_attempts=0,
                  native_residual_kn=float(np.linalg.norm(native_r)),independent_residual_kn=float(np.linalg.norm(residual)),
                  absolute_residual_difference_kn=float(np.linalg.norm(residual-native_r)),
                  independent_residual_over_load=float(np.linalg.norm(residual)/np.linalg.norm(external)),
                  reconstructed_load_difference_kn=float(np.linalg.norm(external-state["applied_nodal_load_kn"])),
                  tangent_relative_frobenius=float(relative_delta),tangent_max_abs_kn_per_m=float(np.max(np.abs(matrix-native_k))),
                  elastic_predictor_relative_to_native=float(np.linalg.norm(elastic_k[np.ix_(global_ids,global_ids)]-native_k)/np.linalg.norm(native_k)),
                  local_stencil_width=h,points_with_mixed_elastic_plastic_probes=sum(min(v["plastic"])!=max(v["plastic"]) for v in stencils.values()),
                  independent_condition_number=float(singular[0]/singular[-1]),
                  independent_min_singular=float(singular[-1]),
                  peak_python_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  peak_native_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                  seconds=time.monotonic()-start)
    print(json.dumps(report,indent=2,allow_nan=False))
    np.savez_compressed(args.work/"paired-operators.npz",native_residual=native_r,independent_residual=residual,
                        native_k=native_k,independent_k=matrix,global_ids=global_ids)

if __name__ == "__main__":
    main()
