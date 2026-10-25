"""Independent analytic overburden/pressure checks of the two-layer example."""
import argparse
import hashlib
import json
from pathlib import Path
import meshio
import numpy as np


def verify(run,solver=None):
    mesh = json.loads((run/"mesh.json").read_text(encoding="utf-8"))
    if mesh.get('tetrahedra10') or any('mohr_coulomb' in m or 'invariant_double_hardening' in m for m in json.loads((run/'ground.json').read_text())['materials']):
        from verify_high_order import verify_ground
        return verify_ground(run,solver=solver)
    state = json.loads((run/"initial-state.json").read_text(encoding="utf-8"))
    config = json.loads((run/"ground.json").read_text(encoding="utf-8"))
    manifest = json.loads((run/"bodies.json").read_text(encoding="utf-8"))
    provenance = json.loads((run/"mesh.provenance.json").read_text(encoding="utf-8"))
    digest = hashlib.sha256((run/"ground.json").read_bytes()).hexdigest()
    assert digest == manifest["ground_sha256"] == provenance["ground_sha256"]
    assert state["ground"] == config
    nodes = np.asarray(mesh["nodes"])
    tetra = np.asarray(mesh["tetrahedra"])
    points = nodes[tetra]
    volumes = np.linalg.det((points[:,1:]-points[:,:1]).transpose(0,2,1))/6
    assert np.all(volumes>0)
    centroids = points.mean(axis=1)
    materials = {m["id"]:m for m in config["materials"]}
    layers = {l["body"]:l for l in config["layers"]}
    region_by_id = {r["id"]:r for r in mesh["regions"]}
    total_weight = 0.
    for region in mesh["regions"]:
        mask = np.asarray(mesh["region_ids"]) == region["id"]
        volume = float(volumes[mask].sum())
        assert np.isclose(volume,region["volume_m3"],rtol=1e-12)
        assert np.all(np.asarray(mesh["material_ids"])[mask] == region["material_id"])
        total_weight += volume*materials[region["material_id"]]["unit_weight_kn_m3"]
    max_stress_error = 0.
    for i,point in enumerate(state["point_states"]):
        assert point["element_index"] == i
        region = region_by_id[mesh["region_ids"][i]]
        layer = layers[region["body"]]
        assert point["material_id"] == mesh["material_ids"][i] == layer["material_id"]
        z = centroids[i,2]
        vertical = config["initialization"]["surface_surcharge_kpa"]
        for r in mesh["regions"]:
            thickness = max(0.,r["bounds_m"][1][2]-max(z,r["bounds_m"][0][2]))
            vertical += thickness*materials[r["material_id"]]["unit_weight_kn_m3"]
        h = config["initialization"]["water_level_m"]
        pore = 0. if h is None else config["initialization"]["water_unit_weight_kn_m3"]*max(h-z,0.)
        stress = np.asarray(point["mcc_state"]["effective_stress_compression_positive_kpa"])
        expected = np.array([layer["k0_x"]*(vertical-pore),layer["k0_y"]*(vertical-pore),vertical-pore,0.,0.,0.])
        max_stress_error = max(max_stress_error,float(np.max(np.abs(stress-expected))))
        assert np.allclose(stress,expected,rtol=1e-12,atol=1e-10)
        assert np.isclose(point["total_stress_compression_positive_kpa"][2],vertical,rtol=1e-12)
        assert np.isclose(point["pore_pressure_kpa"],pore,rtol=1e-12)
        p = stress[:3].mean()
        q2 = 1.5*(np.square(stress[:3]-p).sum()+2*np.square(stress[3:]).sum())
        params = materials[layer["material_id"]]["mcc"]
        pc = point["mcc_state"]["preconsolidation_pressure_kpa"]
        assert q2+params["critical_state_m"]**2*p*(p-pc) <= 1e-8*pc**2
        e = params["specific_volume_ncl_at_reference"]-params["lambda"]*np.log(pc/params["reference_pressure_kpa"])+params["kappa"]*np.log(pc/p)-1
        assert np.isclose(e,point["mcc_state"]["void_ratio"],rtol=1e-12)
        assert np.isclose(e/(1+e),point["porosity"],rtol=1e-12)
    assert max_stress_error<1e-9
    faces = {}
    for i,tet in enumerate(tetra):
        for face in ((0,1,2),(0,1,3),(0,2,3),(1,2,3)):
            key = tuple(sorted(int(tet[a]) for a in face))
            faces.setdefault(key,[]).append(mesh["region_ids"][i])
    interface_count = sum(len(owners)==2 and owners[0]!=owners[1] for owners in faces.values())
    assert interface_count == state["shared_interface_triangles"] > 0
    assert state["true_relative_free_force_residual"]<1e-12
    assert np.linalg.norm(state["equilibrium_error_kn"])/total_weight<1e-12
    vtk = meshio.read(run/"initial-state.vtu")
    assert np.array_equal(vtk.cells_dict["tetra"],tetra)
    assert np.array_equal(vtk.cell_data_dict["material_id"]["tetra"],mesh["material_ids"])
    assert np.allclose(vtk.cell_data_dict["void_ratio"]["tetra"],[p["mcc_state"]["void_ratio"] for p in state["point_states"]])
    report = {"nodes":len(nodes),"tetrahedra":len(tetra),"shared_interface_triangles":interface_count,
              "body_volumes_m3":state["body_volumes_m3"],"total_weight_kn":total_weight,
              "maximum_analytic_stress_error_kpa":max_stress_error,
              "true_relative_free_force_residual":state["true_relative_free_force_residual"],
              "relative_global_force_balance_error":float(np.linalg.norm(state["equilibrium_error_kn"])/total_weight),
              "vtu_roundtrip":"passed"}
    (run/"ground-verification.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    files = sorted(p for p in run.iterdir() if p.is_file() and p.name != "accepted-run.json")
    root = Path(__file__).resolve().parents[1]
    code_paths = [(solver or root/"target/release/soil-fem-3d.exe"),root/"Cargo.lock",root/"rust-toolchain.toml"]
    package = {"schema":"soil-fem-accepted-ground-run/1","files_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
               "implementation_sha256":{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in code_paths}}
    (run/"accepted-run.json").write_text(json.dumps(package,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run",type=Path)
    parser.add_argument("--solver",type=Path)
    args=parser.parse_args();verify(args.run,solver=args.solver)
