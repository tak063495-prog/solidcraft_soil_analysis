"""Independent 3D MCC elastic-branch rate FEM for a CAD Tet10 local cut.

Analytic K0, K=vp/kappa, explicit midpoint integration and dense NumPy solve;
no Rust tangent, return mapping or displacement is used to form the reference.
Uses the same spatial quadrature for material-point comparison. The separate
mechanics verifier compares different volume/surface quadrature rules.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import fem_reference as ref
from verify_high_order import physical,constraints,prepare,verify_stages
from export_stages import export


def reference(mesh,ground,count):
    xyz=np.asarray(mesh['nodes']);cells=ref.cells(mesh);active=np.ones(len(cells),dtype=bool)
    selected=mesh['excavation_sets']['Pit'];active[selected]=False
    used,fixed=constraints(mesh,active,'side_rollers');free=np.flatnonzero(~fixed.reshape(-1))
    if len(free)>2500:raise ValueError('independent dense local reference limited to 2500 free DOFs')
    allq,_=prepare(mesh);q=[allq[ei] for ei in np.flatnonzero(active)]
    b=np.asarray([v['b'] for v in q]);weights=np.asarray([v['weight'] for v in q])
    dofs=(cells[active,:,None]*3+np.arange(3)).reshape(len(q),-1)
    lookup=np.full(xyz.size,-1);lookup[free]=np.arange(len(free));reduced=lookup[dofs];width=dofs.shape[1]
    rows=np.broadcast_to(reduced[:,:,None],(len(q),width,width)).reshape(-1)
    cols=np.broadcast_to(reduced[:,None,:],(len(q),width,width)).reshape(-1)
    mask=(rows>=0)&(cols>=0);rows=rows[mask];cols=cols[mask]
    stresses=[];volumes=[];pcs=[];parameters=[];pores=[]
    for ei,qp in zip(np.flatnonzero(active),q):
        s,pore,material,layer=ref.initial_fields(mesh,ground,ei,qp['position']);m=material['mcc']
        assert layer['consolidation']['definition']=='mcc_pc_over_mean_effective'
        pc=s[:,:3].mean(axis=1)*layer['consolidation']['ratio'];p=s[:,:3].mean(axis=1)
        v=m['specific_volume_ncl_at_reference']-m['lambda']*np.log(pc/m['reference_pressure_kpa'])+m['kappa']*np.log(pc/p)
        stresses.append(s);volumes.append(v);pcs.append(pc);parameters.append(m);pores.append(pore)
    s=np.asarray(stresses);v=np.asarray(volumes);pc=np.asarray(pcs);pore=np.asarray(pores)
    kappa=np.asarray([m['kappa'] for m in parameters])[:,None]
    nu=np.asarray([m['poisson'] for m in parameters])[:,None];m2=np.asarray([m['critical_state_m']**2 for m in parameters])[:,None]
    def internal(stress):
        total=stress.copy();total[:,:,:3]+=pore[:,:,None];out=np.zeros(xyz.size)
        np.add.at(out,dofs.reshape(-1),(-np.einsum('eqai,eqa,eq->ei',b,total,weights,optimize=True)).reshape(-1))
        return out
    pressure=ground['initialization']['surface_surcharge_kpa']
    f,_=physical(mesh,ground,active,pressure,True);f=f.reshape(-1)
    holding=internal(s)-f;holding.reshape(-1,3)[np.setdiff1d(np.arange(len(xyz)),used)]=0
    increment=-holding[free]/count;u=np.zeros(xyz.size);maximum=-1.
    def rate(stress,specific):
        p=stress[:,:,:3].mean(axis=2);dev=stress.copy();dev[:,:,:3]-=p[:,:,None]
        q2=1.5*(np.square(dev[:,:,:3]).sum(axis=2)+2*np.square(dev[:,:,3:]).sum(axis=2))
        yf=(q2+m2*p*(p-pc))/(m2*pc**2)
        assert p.min()>0 and yf.max()<-1e-5,('reference left elastic branch',float(p.min()),float(yf.max()))
        k=specific*p/kappa;g=3*k*(1-2*nu)/(2*(1+nu));lam=k-2*g/3
        d=np.zeros((*p.shape,6,6));d[:,:,:3,:3]=lam[:,:,None,None]
        for a in range(3):d[:,:,a,a]+=2*g;d[:,:,a+3,a+3]=g
        ke=np.einsum('eqai,eqab,eqbj,eq->eij',b,d,b,weights,optimize=True)
        kg=np.zeros((len(free),len(free)));np.add.at(kg,(rows,cols),ke.reshape(-1)[mask])
        du=np.zeros(xyz.size);du[free]=np.linalg.solve(kg,increment)
        de=-np.einsum('eqai,ei->eqa',b,du[dofs],optimize=True)
        ds=np.einsum('eqab,eqb->eqa',d,de,optimize=True)
        return du,ds,-specific*de[:,:,:3].sum(axis=2),float(yf.max())
    for _ in range(count):
        _,ds,dv,y=rate(s,v);du,ds,dv,y2=rate(s+ds/2,v+dv/2)
        u+=du;s+=ds;v+=dv;maximum=max(maximum,y,y2)
    residual=np.linalg.norm((f-internal(s))[free])/max(np.linalg.norm(f),1)
    return {'u':u.reshape(-1,3),'stress':s.reshape(-1,6),'void_ratio':v.reshape(-1)-1,
        'pc':pc.reshape(-1),'maximum_yield':maximum,'force_residual':float(residual),'holding':holding,'free':free,'used':used}


def verify(ground_run,solver,output,steps):
    output.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1]
    mesh_path=ground_run/'mesh.json';mesh=json.loads(mesh_path.read_text());ground=json.loads((ground_run/'ground.json').read_text())
    # This comparison is deliberately well inside MCC's elastic surface. Plastic
    # hardening, cutback and restart have separate theory and state-replay tests.
    for layer in ground['layers']:layer['consolidation']['ratio']=10.
    ground_path=output/'ground.json';ground_path.write_text(json.dumps(ground,indent=2))
    case=json.loads((root/'examples/stages-local-sparse.json').read_text())
    case['settings']['material_maximum_substep']=1e-5
    case['verification_reject_once_stages']=['excavate_pit']
    case_path=output/'stages.json';case_path.write_text(json.dumps(case,indent=2));result_path=output/'stage-results.json'
    subprocess.run([str(solver),'solve-stages',str(mesh_path),str(ground_path),str(case_path),str(result_path)],check=True)
    run=json.loads(result_path.read_text());assert run['complete'] and len(run['frames'])==2
    export(mesh_path,result_path,output/'results');verify_stages(mesh_path,ground_path,result_path,output/'results')
    coarse=reference(mesh,ground,steps);fine=reference(mesh,ground,2*steps)
    state=run['last_committed'];active=np.asarray(state['active_elements']);q,offset=prepare(mesh)
    states=[p for ei in np.flatnonzero(active) for p in state['material_states'][offset[ei]:offset[ei+1]]]
    actual=np.asarray(state['displacement_m']);scale=np.linalg.norm(fine['u'][fine['used']])
    error=np.linalg.norm(actual[fine['used']]-fine['u'][fine['used']])/scale
    refinement=np.linalg.norm(coarse['u'][fine['used']]-fine['u'][fine['used']])/scale
    stress=np.asarray([p['effective_stress_compression_positive_kpa'] for p in states]);e=np.asarray([p['void_ratio'] for p in states])
    stress_error=np.linalg.norm(stress-fine['stress'])/np.linalg.norm(fine['stress'])
    void_error=np.linalg.norm(e-fine['void_ratio'])/np.linalg.norm(fine['void_ratio'])
    assert error<.005 and refinement<1e-4 and stress_error<.001 and void_error<.001,(error,refinement,stress_error,void_error)
    assert fine['force_residual']<1e-6
    assert np.allclose(fine['holding'][fine['free']],np.asarray(run['excavation_audits'][0]['initial_holding_load_kn'])[fine['free']],rtol=1e-9,atol=1e-7)
    assert np.allclose([p['preconsolidation_pressure_kpa'] for p in states],fine['pc'],rtol=1e-12)
    monitors={}
    for name,set_name in case['monitor_sets'].items():
        nodes=np.intersect1d(mesh['node_sets'][set_name],fine['used'])
        monitors[name]={'fem_mean_u_z_m':float(actual[nodes,2].mean()),'reference_mean_u_z_m':float(fine['u'][nodes,2].mean()),
            'fem_minimum_u_z_m':float(actual[nodes,2].min()),'reference_minimum_u_z_m':float(fine['u'][nodes,2].min())}
    assert monitors['excavation_floor']['fem_mean_u_z_m']>0 and monitors['surrounding_surface']['fem_minimum_u_z_m']<0
    assert all(a['committed_state_unchanged_on_rejection'] for a in run['attempts'] if not a['accepted'])
    report={'schema':'soil-fem-independent-quadratic-local/1','reference':'analytic K0 and MCC elastic rates, explicit midpoint, NumPy dense solve; same spatial quadrature',
        'all_active_nodes_three_components_relative_error':float(error),'all_active_material_point_stress_relative_error':float(stress_error),
        'all_active_void_ratio_relative_error':float(void_error),'reference_steps':[steps,2*steps],
        'reference_refinement_error':float(refinement),'independent_force_residual':fine['force_residual'],
        'maximum_normalized_yield':fine['maximum_yield'],'monitors':monitors,
        'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),
        'input_and_reference_sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (mesh_path,ground_path,case_path,Path(__file__),root/'tools/fem_reference.py')}}
    (output/'independent-local-verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--ground-run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--solver',type=Path,default=Path('target/release/soil-fem-3d.exe'));p.add_argument('--steps',type=int,default=100)
    a=p.parse_args();assert a.steps>=20;verify(a.ground_run.resolve(),a.solver.resolve(),a.output_dir.resolve(),a.steps)
