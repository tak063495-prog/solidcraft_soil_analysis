"""Independent 1D rate-equation comparison of the 3D staged MCC FEM.

No Rust return mapping or global stiffness assembly is used in the reference.
Horizontal-plane constraints reduce the verification specimen to an oedometer.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np


def column_mesh():
    nodes = [[float(x),float(y),float(z)] for z in range(11) for y in (0,10) for x in (0,10)]
    tets,materials,regions = [],[],[]
    for z in range(10):
        c = [4*z,4*z+1,4*z+3,4*z+2,4*z+4,4*z+5,4*z+7,4*z+6]
        for t in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)):
            tets.append([c[i] for i in t]);materials.append(1 if z<5 else 2);regions.append(1 if z<5 else 2)
    return {"schema":"soil-fem-mesh/2","units":"m-kN-kPa","axis":"z_up","crs":"LOCAL","origin_m":[0.,0.,0.],
            "nodes":nodes,"tetrahedra":tets,"material_ids":materials,"region_ids":regions,"node_sets":{},
            "regions":[{"id":i+1,"body":name,"material_id":i+1,"bounds_m":[[0.,0.,5.*i],[10.,10.,5.*(i+1)]],
                        "volume_m3":500.,"source_sha256":"a"*64} for i,name in enumerate(("LowerClay","UpperClay"))]}


def reference(initial,ground,count,case):
    states = initial["material_states"]
    stress = np.array([s["effective_stress_compression_positive_kpa"][:3] for s in states])
    pc = np.array([s["preconsolidation_pressure_kpa"] for s in states])
    v = np.array([1+s["void_ratio"] for s in states])
    total_strain = np.zeros(10)
    mats = [ground["materials"][0]["mcc"]]*30+[ground["materials"][1]["mcc"]]*30
    m2=np.array([m["critical_state_m"]**2 for m in mats]);kappa=np.array([m["kappa"] for m in mats])
    hard=np.array([m["lambda"]-m["kappa"] for m in mats]);nu=np.array([m["poisson"] for m in mats])
    water = ground["initialization"]["water_unit_weight_kn_m3"]*5 if ground["initialization"]["water_level_m"] is not None else 0.
    weight = ground["materials"][1]["unit_weight_kn_m3"]*5
    stages=[];pressure=ground['initialization']['surface_surcharge_kpa']
    for stage in case['stages']:
        target=stage.get('target_pressure_kpa',stage.get('remaining_surface_pressure_kpa'))
        stages.append((target-pressure-(weight-water if stage['kind']=='excavate' else 0),5 if stage['kind']=='excavate' else 10))
        pressure=target
    outputs=[]
    for change,active_bands in stages:
        n=6*active_bands
        for _ in range(count):
            s=stress[:n];p=s.mean(axis=1);dev=s-p[:,None];q2=1.5*np.square(dev).sum(axis=1)
            k=v[:n]*p/kappa[:n];g=3*k*(1-2*nu[:n])/(2*(1+nu[:n]));lam=k-2*g/3
            grad=3*dev+m2[:n,None]*(2*p-pc[:n])[:,None]/3
            dg=lam[:,None]*grad.sum(axis=1)[:,None]+2*g[:,None]*grad
            dcol=np.stack((lam,lam,lam+2*g),axis=1)
            numerator=(grad*dcol).sum(axis=1)
            hmod=m2[:n]*p*pc[:n]*v[:n]/hard[:n]*grad.sum(axis=1)
            denominator=(grad*dg).sum(axis=1)+hmod
            yielding=(q2+m2[:n]*p*(p-pc[:n])>=-1e-9*pc[:n]**2)&(numerator*change>0)
            c=dcol[:,2].copy()
            c[yielding]-=dg[yielding,2]*numerator[yielding]/denominator[yielding]
            assert np.all(c>0),"reference stress-controlled path loses positive axial stiffness"
            deps_band=(change/count)/c.reshape(active_bands,6).mean(axis=1)
            deps=np.repeat(deps_band,6)
            multiplier=np.where(yielding,numerator*deps/denominator,0.)
            stress[:n] += dcol*deps[:,None]-dg*multiplier[:,None]
            pc[:n] += pc[:n]*v[:n]/hard[:n]*multiplier*grad.sum(axis=1)
            v[:n] *= np.exp(-deps)
            total_strain[:active_bands] += deps_band
        profile=np.concatenate(([0.],-np.cumsum(total_strain[:active_bands])))
        outputs.append({"profile":profile.copy(),"stress":stress.copy(),"pc":pc.copy(),"void_ratio":v.copy()-1})
    return outputs


def verify(solver,output,small_strain=False):
    output.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[1]
    mesh=column_mesh();ground=json.loads((root/('examples/ground-staged-small.json' if small_strain else 'examples/ground.json')).read_text())
    case=json.loads((root/('examples/stages-small.json' if small_strain else 'examples/stages.json')).read_text())
    case["settings"]["kinematics"]="oedometer"
    case["settings"]["material_maximum_substep"]=1e-5
    case["settings"]["initial_increment"]=0.0625
    for name,data in (("mesh.json",mesh),("ground.json",ground),("case.json",case)):
        (output/name).write_text(json.dumps(data,indent=2),encoding="utf-8")
    subprocess.run([str(solver),"solve-stages",str(output/"mesh.json"),str(output/"ground.json"),str(output/"case.json"),str(output/"result.json")],check=True)
    result=json.loads((output/"result.json").read_text());assert result["complete"]
    coarser=json.loads(json.dumps(case));coarser["settings"]["material_maximum_substep"]=1e-4;coarser["settings"]["initial_increment"]=.125
    (output/"coarse-case.json").write_text(json.dumps(coarser,indent=2),encoding="utf-8")
    subprocess.run([str(solver),"solve-stages",str(output/"mesh.json"),str(output/"ground.json"),str(output/"coarse-case.json"),str(output/"coarse-result.json")],check=True)
    coarse_result=json.loads((output/"coarse-result.json").read_text())
    coarse=reference(result["frames"][0]["state"],ground,5000,case)
    fine=reference(result["frames"][0]["state"],ground,10000,case)
    reports=[]
    for stage,(a,b,frame) in enumerate(zip(coarse,fine,result["frames"][1:])):
        n=len(b["profile"])-1
        u=np.array(frame["state"]["displacement_m"])
        actual=np.array([u[4*z:4*z+4,2].mean() for z in range(n+1)])
        scale=max(np.linalg.norm(b["profile"]),1e-8)
        error=np.linalg.norm(actual-b["profile"])/scale
        uc=np.asarray(coarse_result["frames"][stage+1]["state"]["displacement_m"])
        coarse_profile=np.array([uc[4*z:4*z+4,2].mean() for z in range(n+1)])
        coarse_fem_error=np.linalg.norm(coarse_profile-b["profile"])/scale
        reference_error=np.linalg.norm(a["profile"]-b["profile"])/scale
        assert error<0.005,(frame["name"],error)
        assert error<coarse_fem_error or (error<1e-5 and coarse_fem_error<1e-5)
        assert reference_error<0.001
        states=frame["state"]["material_states"][:6*n]
        actual_pc=np.array([s["preconsolidation_pressure_kpa"] for s in states])
        pc_error=np.linalg.norm(actual_pc-b["pc"][:6*n])/np.linalg.norm(b["pc"][:6*n])
        assert pc_error<0.005
        reports.append({"stage":frame["name"],"surface_displacement_m":actual[-1],
                        "independent_surface_displacement_m":b["profile"][-1],"relative_profile_error":error,
                        "coarse_fem_profile_error":coarse_fem_error,
                        "reference_5000_to_10000_refinement_error":reference_error,"relative_pc_error":pc_error,
                        "true_force_residual":frame["relative_force_residual"]})
    assert reports[0]["surface_displacement_m"]<0
    heave=result["frames"][3]["monitors"]["excavation_floor"]["stage_increment_z_m"]
    independent_heave=fine[2]["profile"][-1]-fine[1]["profile"][5]
    assert heave>0 and independent_heave>0
    assert abs(heave/independent_heave-1)<0.005
    report={"stages":reports,"excavation_heave_increment_m":heave,"independent_heave_increment_m":independent_heave,
            "rollback_attempts":sum(not a["accepted"] for a in result["attempts"])}
    failed_case=json.loads(json.dumps(case));failed_case['stages']=[failed_case['stages'][2]]
    failed_case['settings'].update(maximum_newton_iterations=1,initial_increment=.25,minimum_increment=.25)
    (output/'failed-case.json').write_text(json.dumps(failed_case,indent=2),encoding='utf-8')
    failed_process=subprocess.run([str(solver),'solve-stages',str(output/'mesh.json'),str(output/'ground.json'),str(output/'failed-case.json'),str(output/'failed-result.json')],capture_output=True,text=True)
    assert failed_process.returncode==1
    failed=json.loads((output/'failed-result.json').read_text())
    assert not failed['complete'] and failed['last_committed']==failed['frames'][0]['state']
    from export_stages import export
    try:
        export(output/'mesh.json',output/'failed-result.json',output/'failed-export')
        raise AssertionError('an incomplete calculation was exported as accepted')
    except ValueError as error:
        assert 'completed' in str(error)
    assert not (output/'failed-export').exists()
    report['failure_cli']={'exit_code':failed_process.returncode,'restored_exactly':True,'accepted_export_created':False}
    (output/"verification.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--solver",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument('--small-strain',action='store_true')
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve(),a.small_strain)
