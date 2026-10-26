"""Independent elastic 3D soil + frame solve with all rotations retained.

Frame matrices come from section flexibility/virtual work. Soil uses separate
Gauss-Legendre/Duffy quadrature and SVD constraints. Each installation reference
is formed from the preceding independent solution, never from Rust history.
"""
import argparse
import hashlib
import json
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import fem_reference as ref
from beam_reference import System
from verify_high_order import physical


def verify(run,output):
    if output.exists():raise ValueError('reference output exists')
    mp=run/'ground/mesh.json';gp=run/'ground/ground.json';rp=run/'stage-results.json'
    mesh=json.loads(mp.read_text());ground=ref.normalize_ground(json.loads(gp.read_text()));native=json.loads(rp.read_text());case=native['case']
    plate_case=bool((case.get('beam_system') or {}).get('plates'))
    assert native['complete'] and native['schema']==('soil-fem-stage-result/14' if plate_case else 'soil-fem-stage-result/7')
    watched=[mp,gp,rp,Path(__file__),Path(ref.__file__),Path(__file__).with_name('beam_reference.py')]
    if plate_case:
        solver=run/'solver/soil-fem-3d.exe';binding=run/'solver-binding.json';stages=run/'stages.json'
        digest=hashlib.sha256(solver.read_bytes()).hexdigest()
        assert json.loads(binding.read_text())['solver_sha256']==native['provenance']['solver_binary_sha256']==digest
        assert ref.normalize_stage_case(json.loads(stages.read_text()))==case
        assert native['ground']==ground
        watched.extend([Path(__file__).with_name('plate_reference.py'),solver,binding,stages])
    fingerprints={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in watched}
    assert all(s['kind'] in ('install_beams','remove_beams','boundary','excavate') for s in case['stages'])
    assert not case.get('initially_inactive_bodies') and ground['initialization']['method']=='equilibrate'
    cfg=ground['initialization']['geostatic'];assert cfg['seed']['kind']=='body_states' and cfg['initial_boundary_factor']==0
    for c in case['boundary']['constraints']:assert c['value_m']==0 and not np.asarray(c['gradient']).any()
    xyz=np.asarray(mesh['nodes']);ndof=xyz.size;active=np.ones(len(mesh['tetrahedra']),bool)
    if ndof>6500:raise ValueError('independent direct reference is limited to 6500 soil DOFs')
    materials={m['id']:m for m in ground['materials']};regions={r['id']:r for r in mesh['regions']};seeds=cfg['seed']['states']
    k=np.zeros((ndof,ndof));seed_force=np.zeros(ndof);operators=[];element_operators=[]
    for ei in range(len(active)):
        m=materials[mesh['material_ids'][ei]]['mohr_coulomb'];d=ref.elastic_tensor(m['young_kpa'],m['poisson'])
        seed=seeds[regions[mesh['region_ids'][ei]]['body']];s0=np.asarray(seed['effective_stress_compression_positive_kpa'])
        q=ref.volume(mesh,ei,independent=True);dofs=(q['nodes'][:,None]*3+np.arange(3)).reshape(-1)
        local_k=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        k[np.ix_(dofs,dofs)]+=local_k
        total=np.broadcast_to(s0,(len(q['weight']),6)).copy();h=ground['initialization']['water_level_m']
        if h is not None:total[:,:3]+=np.maximum(h-q['position'][:,2],0)[:,None]*ground['initialization']['water_unit_weight_kn_m3']
        local_seed=-np.einsum('qai,qa,q->i',q['b'],total,q['weight'],optimize=True)
        seed_force[dofs]+=local_seed
        element_operators.append((dofs,local_k,local_seed,float(q['weight'].sum())))
        operators.append((ref.volume(mesh,ei),d,s0,seed['void_ratio']))
    basis,error=ref.staged_constraint_space(mesh,case,{'active_elements':active,'boundary_factor':0,'displacement_m':np.zeros_like(xyz)})
    assert error==0;count=sum(b.shape[1] for b in basis.values());t=np.zeros((ndof,count));start=0
    for n,b in basis.items():t[3*n:3*n+3,start:start+b.shape[1]]=b;start+=b.shape[1]
    beams=System(mesh,case['beam_system']);state=beams.initial_state();u=np.zeros_like(xyz);factor=0.;reports=[]
    pressure=ground['initialization']['surface_surcharge_kpa'];release_reports=[]
    offsets=np.cumsum([0]+[len(qp['weight']) for qp,*_ in operators])
    stress=np.concatenate([np.broadcast_to(s0,(len(qp['weight']),6)).copy() for qp,_,s0,_ in operators])
    void=np.concatenate([np.full(len(qp['weight']),e0) for qp,_,_,e0 in operators])
    for fi,frame in enumerate(native['frames']):
        removed=[]
        if fi:
            stage=case['stages'][fi-1]
            if stage['kind']=='boundary':factor=stage['target_factor']
            elif stage['kind']=='excavate':
                selected=set()
                for body in stage['bodies']:
                    ids=mesh.get('excavation_sets',{}).get(body)
                    if ids is None:
                        region=[r['id'] for r in mesh['regions'] if r['body']==body]
                        assert len(region)==1,'unknown excavation body'
                        ids=[i for i,r in enumerate(mesh['region_ids']) if r==region[0]]
                    assert ids and all(active[i] for i in ids),'empty or previously removed excavation'
                    selected.update(ids)
                removed=sorted(selected);active[removed]=False;assert active.any()
                pressure=stage['remaining_surface_pressure_kpa']
                k=np.zeros((ndof,ndof));seed_force=np.zeros(ndof)
                for ei in np.flatnonzero(active):
                    dofs,local_k,local_seed,_=element_operators[ei]
                    k[np.ix_(dofs,dofs)]+=local_k;seed_force[dofs]+=local_seed
            else:state=beams.change(state,u,stage)
        if removed:
            basis,error=ref.staged_constraint_space(mesh,case,{'active_elements':active,'boundary_factor':factor,'displacement_m':u})
            assert error<1e-10
            count=sum(b.shape[1] for b in basis.values());t=np.zeros((ndof,count));start=0
            for n,b in basis.items():t[3*n:3*n+3,start:start+b.shape[1]]=b;start+=b.shape[1]
            used=set(map(int,np.unique(ref.cells(mesh)[active])))
            assert all(beams.attach[n] is None or beams.attach[n] in used for n in np.flatnonzero(beams.active_nodes(state))), 'active member lost its retained soil attachment'
        kb,rest,fb,eff,extf,extm,info=beams.assemble(state,factor)
        offset,binding,qb,normals=beams.space(state,u,soil_basis=t)
        # Retain beam translations, rotations and released end rotations in
        # the same direct system; no beam Schur complement is used here.
        kk=np.block([[t.T@k@t+binding.T@kb@binding,binding.T@kb@qb],
                     [qb.T@kb@binding,qb.T@kb@qb]])
        fsoil,_=physical(mesh,ground,active,pressure,case.get('pressure_on_original_surface_only',False))
        fsoil+=factor*ref.conditions_load(mesh,case['boundary'],active)
        if removed:
            removed_internal=np.zeros(ndof)
            for ei in removed:
                dofs,local_k,local_seed,_=element_operators[ei]
                removed_internal[dofs]+=local_k@u.reshape(-1)[dofs]+local_seed
            holding=previous_soil_external.reshape(-1)-fsoil.reshape(-1)-removed_internal
            inactive_nodes=np.ones(len(xyz),bool);inactive_nodes[list(used)]=False
            holding.reshape(-1,3)[inactive_nodes]=0.
            actual_audit=next(a for a in native['excavation_audits'] if a['stage']==frame['name'])
            assert np.allclose(holding,actual_audit['initial_holding_load_kn'],atol=1e-7,rtol=1e-8),'independent excavation holding load differs'
            removed_volume=sum(element_operators[ei][3] for ei in removed)
            assert actual_audit['removed_element_count']==len(removed)
            assert np.isclose(actual_audit['removed_volume_m3'],removed_volume,rtol=1e-10,atol=1e-12)
            release_reports.append(dict(stage=frame['name'],removed_elements=removed,removed_volume_m3=removed_volume,
                initial_holding_load_kn=holding.tolist(),maximum_holding_load_error_kn=float(np.max(np.abs(holding-np.asarray(actual_audit['initial_holding_load_kn']))))))
        beam_rhs=fb-rest-kb@offset
        rhs=np.r_[t.T@(fsoil.reshape(-1)-seed_force)+binding.T@beam_rhs,qb.T@beam_rhs]
        solution=np.linalg.solve(kk,rhs);predicted=(t@solution[:count]).reshape(-1,3)
        for n in basis:u[n]=predicted[n]
        beam_d=offset+binding@solution[:count]+qb@solution[count:]
        br=beams.results(state,beam_d,kb,rest,fb,eff,extf,extm,info,normals);state=br['state']
        residual=float(np.linalg.norm(kk@solution-rhs)/max(np.linalg.norm(rhs),1.));assert residual<1e-10
        actual=frame['state'];assert actual['beams']['active_members']==state['active_members']
        assert np.array_equal(actual['active_elements'],active)
        du=float(np.linalg.norm(np.asarray(actual['displacement_m'])-u)/max(np.linalg.norm(u),1e-12));assert du<1e-7,du
        dr=float(np.max(np.abs(np.asarray(actual['beams']['node_rotation_rad'])-state['node_rotation_rad'])))
        db=float(np.max(np.abs(np.asarray(actual['beams']['node_displacement_m'])-state['node_displacement_m'])))
        refs=float(np.max(np.abs(np.asarray(actual['beams']['reference_global_kinematics'])-state['reference_global_kinematics']),initial=0.))
        if plate_case:
            assert actual['beams']['active_plates']==state['active_plates']
            plate_reference_error=float(np.max(np.abs(np.asarray(actual['beams']['plate_reference_global_kinematics'])-state['plate_reference_global_kinematics']),initial=0.))
            refs=max(refs,plate_reference_error)
        assert dr<1e-9 and db<1e-9 and refs<1e-9,(dr,db,refs)
        for ei,(qp,d,s0,e0) in enumerate(operators):
            if not active[ei]:continue
            strain=np.einsum('qai,i->qa',qp['b'],u[qp['nodes']].reshape(-1))
            stress[offsets[ei]:offsets[ei+1]]=s0-strain@d.T
            void[offsets[ei]:offsets[ei+1]]=(1+e0)*np.exp(strain[:,:3].sum(axis=1))-1
        points=actual['material_states'];assert all(p['equivalent_plastic_strain']==0 for p in points),'elastic reference entered MC plasticity'
        stress=np.asarray(stress);ds=float(np.linalg.norm(np.asarray([p['effective_stress_compression_positive_kpa'] for p in points])-stress)/np.linalg.norm(stress))
        de=float(np.max(np.abs(np.asarray([p['void_ratio'] for p in points])-void)));assert ds<1e-7 and de<1e-8,(ds,de)
        balance=(k@u.reshape(-1)+seed_force-fsoil.reshape(-1)).reshape(-1,3)
        for node,value in br['interface_residual'].items():balance[node]+=value
        support=np.zeros_like(xyz)
        for node,b in basis.items():support[node]=balance[node]-b@(b.T@balance[node])
        assert np.allclose(support.sum(axis=0)+br['support_force'],frame['support_reaction_kn'],atol=1e-5,rtol=1e-9)
        moment=np.cross(xyz,fsoil+support).sum(axis=0)+br['external_moment']+br['support_moment']
        assert np.linalg.norm(moment)<1e-6
        error=0.
        for member in frame['beam_results']['members']:
            value=np.asarray(member['local_end_actions']);expected=br['members'][member['id']]['local_end_actions']
            error=max(error,float(np.max(np.abs(value-expected))))
            assert np.allclose(value,expected,atol=1e-6,rtol=1e-7),(value,expected)
        if plate_case:
            assert {p['id'] for p in frame['beam_results'].get('plates',[])}==set(br['plates'])
            for plate in frame['beam_results'].get('plates',[]):
                expected=br['plates'][plate['id']]
                for key in ('area_m2','local_axes','membrane_energy_kn_m','bending_energy_kn_m','drilling_penalty_energy_kn_m'):
                    assert np.allclose(plate[key],expected[key],atol=1e-8,rtol=1e-7),(plate['id'],key)
                for actual_point,reference_point in zip(plate['samples'],expected['samples'],strict=True):
                    for key in reference_point:
                        assert np.allclose(actual_point[key],reference_point[key],atol=1e-8,rtol=1e-7),(plate['id'],key)
        reports.append({'frame':frame['name'],'soil_free_dofs':count,'retained_frame_unknowns':qb.shape[1],
                        'relative_joint_system_residual':residual,'relative_3d_displacement_error':du,
                        'maximum_beam_displacement_error_m':db,'maximum_rotation_error_rad':dr,
                        'maximum_installation_reference_error':refs,'maximum_end_action_error':error,
                        'relative_all_point_stress_error':ds,'maximum_void_ratio_error':de,
                        'global_moment_balance_kn_m':moment.tolist()})
        previous_soil_external=fsoil.copy()
    report={'schema':'soil-fem-independent-beam-coupling/1','reference':'full soil + frame direct NumPy solve, input seed and independently propagated installation references; all rotations and releases retained; section flexibility by virtual work; GL125 soil integration; SVD constraints',
            'frames':reports,'solver_sha256':native['provenance']['solver_binary_sha256'],
            'sha256':fingerprints}
    if plate_case:
        report.update(schema='soil-fem-independent-plate-beam-coupling/1',
            reference='Full soil + plate/frame direct solve; 12 independent polynomial DKT constraints and GL25 surface integration; GL125 soil quadrature; SVD supports; no Rust solved field or installation reference generates the solution',
            plate_samples_and_energies_checked=True,spatial_convergence_accepted=False,goal_complete=False,p6_complete=False)
    if release_reports:
        report.update(excavation_release_audits=release_reports,
            removed_soil_displacement_and_gp_histories_independently_frozen=True,
            scope='Elastic soil + linear plates/frames, supported retained attachments; not plastic excavation or general soil/plate interface convergence')
    assert fingerprints=={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in watched},'reference inputs or sources changed'
    with output.open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();verify(args.run.resolve(),args.output.resolve())
