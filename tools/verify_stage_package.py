"""Independent numpy assembly of accepted material forces, loads and reactions."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def result_schema(path):
    # Header routing must not retain a second copy of every GP history.
    import codecs,ijson
    from streamed_stage_source import CheckedReader
    from ijson.common import JSONError
    with path.open('rb') as handle:
        if handle.read(3)!=codecs.BOM_UTF8:handle.seek(0)
        try:
            for prefix,event,value in ijson.parse(CheckedReader(handle,path.stat().st_size),use_float=True):
                if prefix=='schema' and event=='string':return value
        except JSONError as error:raise ValueError('invalid result header JSON') from error
    raise ValueError('result schema string required')


def verify(mesh_path,ground_path,result_path,export_dir,solver=None):
    mesh=json.loads(mesh_path.read_text(encoding="utf-8"));ground=json.loads(ground_path.read_text(encoding="utf-8"))
    schema=result_schema(result_path)
    if schema in ('soil-fem-mixed-mcc-stage-result/1','soil-fem-mixed-mcc-stage-result/2','soil-fem-projected-mcc-stage-result/1','soil-fem-projected-mcc-stage-result/2'):
        from verify_projected_stage_package import verify as verify_projected
        return verify_projected(mesh_path,ground_path,result_path,export_dir,solver=solver)
    result=json.loads(result_path.read_text(encoding="utf-8"))
    if result.get('schema')!=schema:raise ValueError('result schema changed or duplicated')
    if mesh.get('tetrahedra10') or result.get('case',{}).get('boundary') is not None or result.get('schema') in ('soil-fem-stage-result/3','soil-fem-stage-result/4','soil-fem-stage-result/5','soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/13','soil-fem-stage-result/14'):
        from verify_high_order import verify_stages
        return verify_stages(mesh_path,ground_path,result_path,export_dir,solver=solver)
    assert result["complete"] and result["ground"]==ground
    xyz=np.asarray(mesh["nodes"]);tets=np.asarray(mesh["tetrahedra"]);bounds=[xyz.min(axis=0),xyz.max(axis=0)]
    material={m["id"]:m for m in ground["materials"]}
    stage_case=json.loads((result_path.parent/"stages.json").read_text(encoding="utf-8"))
    from fem_reference import normalize_stage_case
    normalize_stage_case(stage_case)
    assert stage_case==result["case"]
    def pore(z):
        h=ground["initialization"]["water_level_m"]
        return 0. if h is None else max(h-z,0.)*ground["initialization"]["water_unit_weight_kn_m3"]
    reports=[]
    for index,frame in enumerate(result["frames"]):
        state=frame["state"];active=np.asarray(state["active_elements"]);f=np.zeros_like(xyz);internal=np.zeros_like(xyz)
        used=np.unique(tets[active]);faces={};max_yield=-np.inf
        for i in np.flatnonzero(active):
            tet=tets[i];p=xyz[tet];j=(p[1:]-p[0]).T;volume=np.linalg.det(j)/6
            gradients=np.linalg.inv(j);g=np.vstack((-gradients.sum(axis=0),gradients))
            params=material[mesh["material_ids"][i]]
            s=np.asarray(state["material_states"][i]["effective_stress_compression_positive_kpa"])
            mean=s[:3].mean();q2=1.5*(np.square(s[:3]-mean).sum()+2*np.square(s[3:]).sum())
            pc=state["material_states"][i]["preconsolidation_pressure_kpa"];e=state["material_states"][i]["void_ratio"]
            m2=params["mcc"]["critical_state_m"]**2
            yf=(q2+m2*mean*(mean-pc))/(m2*pc**2);max_yield=max(max_yield,yf)
            assert np.isfinite(s).all() and pc>0 and e>0 and mean>0 and yf<=1e-8
            total=s.copy();total[:3]+=pore(p[:,2].mean())
            stress=np.array([[total[0],total[3],total[5]],[total[3],total[1],total[4]],[total[5],total[4],total[2]]])
            for node,gradient in zip(tet,g):
                internal[node]-=volume*(stress@gradient)
                f[node,2]-=volume*params["unit_weight_kn_m3"]/4
            for face in ((0,1,2),(0,1,3),(0,2,3),(1,2,3)):
                key=tuple(sorted(int(tet[a]) for a in face));faces.setdefault(key,[]).append(p.mean(axis=0))
        all_faces={}
        for ei,tet in enumerate(tets):
            for face in ((0,1,2),(0,1,3),(0,2,3),(1,2,3)):
                key=tuple(sorted(int(tet[a]) for a in face));all_faces.setdefault(key,[]).append(ei)
        cut_area=0.;cut_water=np.zeros(3);original_surface_area=0.;cut_faces=0
        for face,owners in faces.items():
            if len(owners)!=1:continue
            a,b,c=xyz[list(face)];normal=np.cross(b-a,c-a);center=(a+b+c)/3
            if normal@(center-owners[0])<0:normal=-normal
            up=normal[2]>0 and abs(normal[0])+abs(normal[1])<1e-8*np.linalg.norm(normal)
            cut=len(all_faces[face])==2
            if cut:cut_area+=np.linalg.norm(normal)/2;cut_faces+=1
            if cut or up:
                wp=np.array([pore(z) for z in xyz[list(face),2]])
                for a,node in enumerate(face):
                    traction=-normal*(wp.sum()+wp[a])/24
                    f[node]+=traction
                    if cut:cut_water+=traction
            if up and (not result['case'].get('pressure_on_original_surface_only',False) or abs(center[2]-bounds[1][2])<1e-8):
                f[list(face),2]-=np.linalg.norm(normal)/2*state['surface_pressure_kpa']/3
                original_surface_area+=np.linalg.norm(normal)/2
        assert np.allclose(f.reshape(-1),state["applied_nodal_load_kn"],rtol=1e-11,atol=1e-8)
        assert np.linalg.norm(state["excavation_release_load_kn"])<1e-9
        fixed=np.ones_like(xyz,dtype=bool);fixed[used]=False
        if result["case"]["settings"]["kinematics"]=="oedometer":
            fixed[used,:2]=True
        else:
            for axis in (0,1):
                fixed[used,axis]=np.isclose(xyz[used,axis],bounds[0][axis],atol=1e-8,rtol=0)|np.isclose(xyz[used,axis],bounds[1][axis],atol=1e-8,rtol=0)
        fixed[used,2]=np.isclose(xyz[used,2],bounds[0][2],atol=1e-8,rtol=0)
        residual=f-internal
        if result["case"]["settings"]["kinematics"]=="oedometer":
            sums={}
            for n in used:
                if not fixed[n,2]:sums[xyz[n,2]]=sums.get(xyz[n,2],0.)+residual[n,2]
            reduced=np.array(list(sums.values()))
        else:reduced=residual[~fixed]
        relative=np.linalg.norm(reduced)/max(np.linalg.norm(f),1.)
        assert relative<=result["case"]["settings"]["relative_force_tolerance"]*1.01
        reaction=np.where(fixed,internal-f,0.);reaction[np.setdiff1d(np.arange(len(xyz)),used)]=0
        applied=f.sum(axis=0);support=reaction.sum(axis=0)
        assert np.allclose(applied,frame["applied_force_kn"],rtol=1e-12,atol=1e-8)
        assert np.allclose(support,frame["support_reaction_kn"],rtol=1e-10,atol=1e-7)
        relative_moment=None
        if result["case"]["settings"]["kinematics"]=="side_rollers":
            moment=np.cross(xyz,f+reaction).sum(axis=0)
            relative_moment=np.linalg.norm(moment)/max(np.linalg.norm(f)*np.linalg.norm(bounds[1]-bounds[0]),1.)
            assert relative_moment<1e-7
        if index==0:assert not np.asarray(state["displacement_m"]).any()
        if index>0:
            previous=result["frames"][index-1]["state"]
            for ei in np.flatnonzero(~active):
                assert state["material_states"][ei]==previous["material_states"][ei]
        reports.append({"stage":frame["name"],"active_elements":int(active.sum()),"independent_force_residual":relative,
                        "force_balance_error_kn":(applied+support).tolist(),"relative_moment_error":relative_moment,
                        "maximum_normalized_yield":max_yield,"maximum_absolute_engineering_strain_from_k0":frame["maximum_absolute_engineering_strain_from_k0"],
                        "cut_surface_area_m2":cut_area,"cut_face_count":cut_faces,"cut_water_force_kn":cut_water.tolist(),"loaded_surface_area_m2":original_surface_area})
    assert all(a["committed_state_unchanged_on_rejection"] for a in result["attempts"] if not a["accepted"])
    export=json.loads((export_dir/"export-verification.json").read_text())
    assert export["result_sha256"]==hashlib.sha256(result_path.read_bytes()).hexdigest()
    assert export["mesh_sha256"]==hashlib.sha256(mesh_path.read_bytes()).hexdigest()
    report={"schema":"soil-fem-stage-verification/1","frames":reports,"vtu_roundtrip":export["vtu_roundtrip"]}
    (result_path.parent/"stage-verification.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    root=Path(__file__).resolve().parents[1]
    producer=result.get('provenance',{}).get('solver_binary_sha256')
    assert producer==hashlib.sha256(((solver or root/'target/release/soil-fem-3d.exe')).read_bytes()).hexdigest(),'result produced by a different binary'
    package={"schema":"soil-fem-accepted-staged-run/1","inputs_and_results_sha256":{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in [mesh_path,ground_path,result_path,result_path.parent/"stages.json"]},
             "solver_binary_sha256":hashlib.sha256((solver or root/"target/release/soil-fem-3d.exe").read_bytes()).hexdigest()}
    (result_path.parent/"accepted-stage-run.json").write_text(json.dumps(package,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mesh",type=Path);p.add_argument("ground",type=Path);p.add_argument("result",type=Path);p.add_argument("exports",type=Path)
    p.add_argument("--solver",type=Path)
    a=p.parse_args();verify(a.mesh,a.ground,a.result,a.exports,solver=a.solver)
