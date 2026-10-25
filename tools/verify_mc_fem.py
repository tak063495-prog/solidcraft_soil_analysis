"""Separate CLI processes: analytic two-material MC column, force audit/restart."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from fem_reference import EDGES
from export_stages import export
from verify_high_order import verify_stages
from verify_checkpoint import invoke


def theory(total):
    nu=.1;phi=np.sin(np.deg2rad(30));psi=np.sin(np.deg2rad(5));a=(1-nu)/((1+nu)*(1-2*nu));b=nu/((1+nu)*(1-2*nu))
    m=np.array([-(1+psi)/2,-(1+psi)/2,1-psi]);n=np.array([-(1+phi)/2,-(1+phi)/2,1-phi]);d=np.full((3,3),b);np.fill_diagonal(d,a)
    dm=d@m;den=n@dm;drive=(1-phi)*a-(1+phi)*b;initial=(1-phi)*100-(1+phi)*50-20*np.cos(np.deg2rad(30))
    def strain(s,e):
        ey=-initial/(drive*e);sy=100+a*e*ey
        return (s-100)/(a*e) if s<=sy else ey+(s-sy)/(e*(a-dm[2]*drive/den))
    lo=100;hi=1000
    for _ in range(100):
        mid=(lo+hi)/2
        if (strain(mid,10000)+strain(mid,20000))/2<total:lo=mid
        else:hi=mid
    stress=(lo+hi)/2
    return stress,np.array([strain(stress,10000),strain(stress,20000)])


def elevate(mesh):
    mesh=json.loads(json.dumps(mesh));cache={};out=[]
    for t in mesh['tetrahedra']:
        mids=[]
        for a,b in EDGES:
            edge=tuple(sorted((t[a],t[b])))
            if edge not in cache:
                cache[edge]=len(mesh['nodes']);mesh['nodes'].append(((np.array(mesh['nodes'][edge[0]])+mesh['nodes'][edge[1]])/2).tolist())
            mids.append(cache[edge])
        out.append(t+mids)
    for name,(axis,value) in {'xmin':(0,0),'xmax':(0,1),'ymin':(1,0),'ymax':(1,1),'bottom':(2,0),'top':(2,1)}.items():
        mesh['node_sets'][name]=[i for i,x in enumerate(mesh['nodes']) if abs(x[axis]-value)<1e-12]
    mesh['node_sets']['all']=list(range(len(mesh['nodes'])))
    mesh.update(schema='soil-fem-mesh/3',tetrahedra10=out,quadrature_order=3)
    return mesh


def verify(solver,output):
    output.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1];reports=[]
    expected,strains=theory(.02);nu=.1;a=(1-nu)/((1+nu)*(1-2*nu));unload_stress=expected-.01/(.5/(a*10000)+.5/(a*20000))
    for quadratic in (False,True):
        folder=output/('tet10' if quadratic else 'tet4');folder.mkdir();mesh=json.loads((root/'examples/mc-column-mesh.json').read_text())
        if quadratic:mesh=elevate(mesh)
        mp=folder/'mesh.json';mp.write_text(json.dumps(mesh));gp=folder/'ground.json';gp.write_bytes((root/'examples/ground-mc-column.json').read_bytes())
        cp=folder/'stages.json';cp.write_bytes((root/'examples/stages-mc-column.json').read_bytes());rp=folder/'stage-results.json'
        inputs=[str(mp),str(gp),str(cp)]
        invoke([str(solver),'solve-stages',*inputs,str(rp)],folder/'solve.log')
        result=json.loads(rp.read_text());assert result['complete'];export(mp,rp,folder/'results');verify_stages(mp,gp,rp,folder/'results',solver=solver)
        assert any(not a['accepted'] and a['rejected_trial_plastic_changed'] and a['committed_state_unchanged_on_rejection'] for a in result['attempts'])
        errors=[]
        for frame,target in zip(result['frames'][1:],(expected,unload_stress)):
            stresses=np.array([p['effective_stress_compression_positive_kpa'] for p in frame['state']['material_states']]);error=float(np.max(np.abs(stresses[:,2]-target)))
            assert error<2e-5,(error,target);errors.append(error)
        z=np.array(mesh['nodes'])[:,2];u=np.array(result['frames'][1]['state']['displacement_m']);expected_u=-np.where(z<=.5,z*strains[0],.5*strains[0]+(z-.5)*strains[1])
        assert np.max(np.abs(u[:,2]-expected_u))<2e-8
        assert np.max(np.abs(u[:,:2]))<2e-8
        before=result['frames'][1]['state']['material_states'];after=result['frames'][2]['state']['material_states']
        assert all(np.allclose(x['plastic_strain_compression_positive'],y['plastic_strain_compression_positive'],rtol=1e-8,atol=1e-9) for x,y in zip(before,after))
        if quadratic:
            check=folder/'checkpoint.json';resumed=folder/'resumed.json'
            invoke([str(solver),'checkpoint-stages',*inputs,'2',str(check)],folder/'checkpoint.log')
            saved=json.loads(check.read_text());assert any(p['equivalent_plastic_strain']>1e-6 for p in saved['run']['last_committed']['material_states'])
            invoke([str(solver),'resume-stages',*inputs,str(check),str(resumed)],folder/'resume.log')
            resumed_result=json.loads(resumed.read_text());assert result['last_committed']==resumed_result['last_committed'] and result['frames']==resumed_result['frames']
        reports.append({'element_order':2 if quadratic else 1,'integration_points':len(result['integration_points']),
            'theoretical_axial_stress_kpa':[expected,unload_stress],'maximum_absolute_stress_error_kpa':errors,
            'maximum_absolute_vertical_displacement_error_m':float(np.max(np.abs(u[:,2]-expected_u))),
            'accepted_increments':sum(a['accepted'] for a in result['attempts']),
            'restored_on_rejection':all(a['committed_state_unchanged_on_rejection'] for a in result['attempts'] if not a['accepted'])})
    report={'schema':'soil-fem-independent-mc-column/1','method':'closed-form uniaxial-strain MC compression-edge compliance and series compatibility, followed by elastic unloading; independent 3D force audit',
        'columns':reports,'separate_process_tet10_plastic_restart_exact':True,'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest()}
    report['artifacts_sha256']={str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in output.rglob('*.json')}
    for p in (Path(__file__),root/'tools/fem_reference.py',root/'tools/verify_high_order.py',root/'tools/mc_reference.py'):report['artifacts_sha256'][str(p.resolve())]=hashlib.sha256(p.read_bytes()).hexdigest()
    (output/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='artifacts_sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,default=Path('target/release/soil-fem-3d.exe'));p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
