"""Independent 3D rate-form FEM on the identical accepted CAD mesh.

Uses analytic MCC elastic-branch rates with explicit midpoint integration;
no Rust tangent, return mapping, or displacement is used by the reference.
This benchmark must stay strictly inside the yield surface. Plastic local-cut
history is separately exercised by tests/local_excavation.rs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import os
# Small dense systems run more predictably without BLAS thread oversubscription.
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from export_stages import export


def operators(mesh, active, ground):
    xyz=np.asarray(mesh['nodes']);tets=np.asarray(mesh['tetrahedra'])[active]
    p=xyz[tets];j=(p[:,1:]-p[:,:1]).transpose(0,2,1)
    volume=np.linalg.det(j)/6;g=np.linalg.inv(j);g=np.concatenate((-g.sum(axis=1)[:,None],g),axis=1)
    b=np.zeros((len(tets),6,12))
    for a in range(4):
        x,y,z=g[:,a].T;i=3*a
        b[:,0,i]=x;b[:,1,i+1]=y;b[:,2,i+2]=z
        b[:,3,i]=y;b[:,3,i+1]=x;b[:,4,i+1]=z;b[:,4,i+2]=y;b[:,5,i]=z;b[:,5,i+2]=x
    used=np.unique(tets);fixed=np.ones((len(xyz),3),dtype=bool);fixed[used]=False
    lo=xyz.min(axis=0);hi=xyz.max(axis=0)
    for a in (0,1):fixed[used,a]=np.isclose(xyz[used,a],lo[a],atol=1e-8,rtol=0)|np.isclose(xyz[used,a],hi[a],atol=1e-8,rtol=0)
    fixed[used,2]=np.isclose(xyz[used,2],lo[2],atol=1e-8,rtol=0)
    free=np.flatnonzero(~fixed.reshape(-1));dofs=(tets[:,:,None]*3+np.arange(3)).reshape(-1,12)
    lookup=np.full(xyz.size,-1);lookup[free]=np.arange(len(free));reduced=lookup[dofs]
    rows=np.broadcast_to(reduced[:,:,None],(len(tets),12,12)).reshape(-1)
    cols=np.broadcast_to(reduced[:,None,:],(len(tets),12,12)).reshape(-1)
    mask=(rows>=0)&(cols>=0)
    mats={m['id']:m['mcc'] for m in ground['materials']}
    parameters=[mats[i] for i in np.asarray(mesh['material_ids'])[active]]
    return xyz,tets,volume,b,free,dofs,rows[mask],cols[mask],mask,parameters


def reference(mesh,ground,initial,final,holding,count):
    active=np.asarray(final['active_elements']);xyz,tets,vol,b,free,dofs,rows,cols,mask,mats=operators(mesh,active,ground)
    # Analytic K0 fields from depth and layer weights, independently of Rust state.
    regions={r['id']:r for r in mesh['regions']};layers={r['body']:r for r in ground['layers']};materials={m['id']:m for m in ground['materials']}
    all_centers=xyz[np.asarray(mesh['tetrahedra'])].mean(axis=1);centers=all_centers[active]
    s=np.zeros((len(tets),6));pc=np.empty(len(tets));v=np.empty(len(tets))
    for i,(ei,z) in enumerate(zip(np.flatnonzero(active),centers[:,2])):
        total=ground['initialization']['surface_surcharge_kpa']
        for r in mesh['regions']:
            total+=max(r['bounds_m'][1][2]-max(z,r['bounds_m'][0][2]),0)*materials[r['material_id']]['unit_weight_kn_m3']
        h=ground['initialization']['water_level_m'];pore=0 if h is None else max(h-z,0)*ground['initialization']['water_unit_weight_kn_m3']
        layer=layers[regions[mesh['region_ids'][ei]]['body']];m=mats[i]
        s[i,:3]=np.array([layer['k0_x'],layer['k0_y'],1])*(total-pore)
        mean=s[i,:3].mean();assert layer['consolidation']['definition']=='mcc_pc_over_mean_effective'
        pc[i]=mean*layer['consolidation']['ratio']
        v[i]=m['specific_volume_ncl_at_reference']-m['lambda']*np.log(pc[i]/m['reference_pressure_kpa'])+m['kappa']*np.log(pc[i]/mean)
    expected=np.asarray([q['effective_stress_compression_positive_kpa'] for q,a in zip(initial['material_states'],active) if a])
    assert np.allclose(s,expected,rtol=1e-12,atol=1e-10)
    kappa=np.array([m['kappa'] for m in mats]);nu=np.array([m['poisson'] for m in mats]);m2=np.array([m['critical_state_m']**2 for m in mats])
    u=np.zeros(xyz.size);increment=-np.asarray(holding)[free]/count
    def rate(stress,specific_volume):
        p=stress[:,:3].mean(axis=1);dev=stress.copy();dev[:,:3]-=p[:,None]
        q2=1.5*(np.square(dev[:,:3]).sum(axis=1)+2*np.square(dev[:,3:]).sum(axis=1))
        yf=(q2+m2*p*(p-pc))/(m2*pc**2)
        assert np.max(yf)<-1e-5,('benchmark left elastic MCC branch',float(np.max(yf)))
        k=specific_volume*p/kappa;g=3*k*(1-2*nu)/(2*(1+nu));lam=k-2*g/3
        d=np.zeros((len(tets),6,6));d[:,:3,:3]=lam[:,None,None]
        for a in range(3):d[:,a,a]+=2*g;d[:,a+3,a+3]=g
        ke=np.einsum('eai,eab,ebj,e->eij',b,d,b,vol,optimize=True)
        global_k=np.zeros((len(free),len(free)));np.add.at(global_k,(rows,cols),ke.reshape(-1)[mask])
        du=np.zeros(xyz.size);du[free]=np.linalg.solve(global_k,increment)
        de=-np.einsum('eai,ei->ea',b,du[dofs]);ds=np.einsum('eab,eb->ea',d,de)
        return du,ds,-specific_volume*de[:,:3].sum(axis=1),float(np.max(yf))
    maximum_yield=-1.
    for _ in range(count):
        _,ds,dv,y=rate(s,v)
        du,ds,dv,y2=rate(s+ds/2,v+dv/2)
        u+=du;s+=ds;v+=dv;maximum_yield=max(maximum_yield,y,y2)
    # Independent integration drift / equilibrium check at the endpoint.
    internal=np.zeros(xyz.size)
    h=ground['initialization']['water_level_m'];total=s.copy()
    if h is not None:total[:,:3]+=np.maximum(h-centers[:,2],0)[:,None]*ground['initialization']['water_unit_weight_kn_m3']
    fe=-vol[:,None]*np.einsum('eai,ea->ei',b,total);np.add.at(internal,dofs.reshape(-1),fe.reshape(-1))
    residual=np.linalg.norm((np.asarray(final['applied_nodal_load_kn'])-internal)[free])/max(np.linalg.norm(final['applied_nodal_load_kn']),1)
    return {'u':u.reshape(-1,3),'stress':s,'void_ratio':v-1,'maximum_yield':maximum_yield,'force_residual':residual}


def verify(run,solver,counts=(100,200)):
    mesh_path=run/'ground/mesh.json';ground_path=run/'ground/ground.json';result_path=run/'stage-results.json'
    mesh=json.loads(mesh_path.read_text());ground=json.loads(ground_path.read_text());result=json.loads(result_path.read_text())
    assert result['complete'] and len(result['frames'])==2
    for name in ('reference-refined-stages.json','reference-refined-result.json','independent-local-verification.json'):
        if (run/name).exists():raise ValueError(f'comparison output exists: {run/name}')
    initial=result['frames'][0]['state'];end=result['frames'][1]['state'];audit=result['excavation_audits'][0]
    # Reassemble lambda=0 holding forces independently: retain old material stress,
    # subtract final physical load. No audit vector is fed into the reference solve.
    active=np.asarray(end['active_elements']);xyz,tets,vol,b,free,dofs,*_=operators(mesh,active,ground)
    old=np.asarray([q['effective_stress_compression_positive_kpa'] for q,a in zip(initial['material_states'],active) if a])
    total=old.copy();h=ground['initialization']['water_level_m']
    if h is not None:total[:,:3]+=np.maximum(h-xyz[tets].mean(axis=1)[:,2],0)[:,None]*ground['initialization']['water_unit_weight_kn_m3']
    internal=np.zeros(xyz.size);np.add.at(internal,dofs.reshape(-1),(-vol[:,None]*np.einsum('eai,ea->ei',b,total)).reshape(-1))
    holding=internal-np.asarray(end['applied_nodal_load_kn']);used=np.unique(tets);orphan=np.setdiff1d(np.arange(len(xyz)),used);holding.reshape(-1,3)[orphan]=0
    assert np.allclose(holding[free],np.asarray(audit['initial_holding_load_kn'])[free],rtol=1e-9,atol=1e-7)
    assert audit['zero_release_relative_force_residual']<1e-8
    coarse=reference(mesh,ground,initial,end,holding,counts[0]);fine=reference(mesh,ground,initial,end,holding,counts[1])
    # Refine the Rust integration independently as well.
    case=json.loads(json.dumps(result['case']));case['settings']['initial_increment']=0.125;case['settings']['material_maximum_substep']=1e-5
    case_path=run/'reference-refined-stages.json';fine_path=run/'reference-refined-result.json'
    case_path.write_text(json.dumps(case,indent=2));subprocess.run([str(solver),'solve-stages',str(mesh_path),str(ground_path),str(case_path),str(fine_path)],check=True)
    refined=json.loads(fine_path.read_text());assert refined['complete'];actual=np.asarray(refined['last_committed']['displacement_m'])
    export(mesh_path,fine_path,run/'reference-refined-exports')
    coarse_fem=np.asarray(end['displacement_m']);scale=np.linalg.norm(fine['u'][used]);error=np.linalg.norm(actual[used]-fine['u'][used])/scale
    coarse_error=np.linalg.norm(coarse_fem[used]-fine['u'][used])/scale
    refinement=np.linalg.norm(fine['u'][used]-coarse['u'][used])/scale
    assert error<0.005,(error,coarse_error);assert error<coarse_error or coarse_error<1e-5
    assert refinement<1e-4 and fine['force_residual']<1e-6
    states=[q for q,a in zip(refined['last_committed']['material_states'],active) if a]
    actual_stress=np.array([q['effective_stress_compression_positive_kpa'] for q in states]);actual_e=np.array([q['void_ratio'] for q in states])
    stress_error=np.linalg.norm(actual_stress-fine['stress'])/np.linalg.norm(fine['stress'])
    void_ratio_error=np.linalg.norm(actual_e-fine['void_ratio'])/np.linalg.norm(fine['void_ratio'])
    assert stress_error<0.001 and void_ratio_error<0.001
    reports={}
    for name,set_name in result['case']['monitor_sets'].items():
        nodes=np.intersect1d(mesh['node_sets'][set_name],used)
        reports[name]={'fem_mean_u_z_m':float(actual[nodes,2].mean()),'independent_mean_u_z_m':float(fine['u'][nodes,2].mean()),'fem_minimum_u_z_m':float(actual[nodes,2].min()),'independent_minimum_u_z_m':float(fine['u'][nodes,2].min())}
    assert reports['excavation_floor']['fem_mean_u_z_m']>0
    assert reports['surrounding_surface']['fem_minimum_u_z_m']<0
    for ei,a in enumerate(active):
        if a:assert refined['last_committed']['material_states'][ei]['preconsolidation_pressure_kpa']==initial['material_states'][ei]['preconsolidation_pressure_kpa']
    report={'schema':'soil-fem-independent-local-3d/1','method':'independent analytic K0 + analytic MCC elastic rate + explicit midpoint global FEM','relative_3d_displacement_error':error,'coarse_fem_error':coarse_error,'relative_stress_error':stress_error,'relative_void_ratio_error':void_ratio_error,'reference_counts':counts,'reference_refinement_error':refinement,'independent_force_residual':fine['force_residual'],'maximum_normalized_yield':fine['maximum_yield'],'monitors':reports,'compared_fields':'all active nodes, all three displacement components; all active element stresses and void ratios; same spatial discretization','input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (mesh_path,ground_path,result_path)}}
    report['comparison_artifacts_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (case_path,fine_path,Path(__file__))}
    (run/'independent-local-verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--solver',type=Path,default=Path('target/release/soil-fem-3d.exe'));p.add_argument('--steps',type=int,default=100)
    a=p.parse_args();verify(a.run,a.solver.resolve(),(a.steps,2*a.steps))
