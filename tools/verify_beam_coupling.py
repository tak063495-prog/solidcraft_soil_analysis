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
    assert native['complete'] and native['schema']=='soil-fem-stage-result/7'
    assert all(s['kind'] in ('install_beams','remove_beams','boundary') for s in case['stages'])
    assert not case.get('initially_inactive_bodies') and ground['initialization']['method']=='equilibrate'
    cfg=ground['initialization']['geostatic'];assert cfg['seed']['kind']=='body_states' and cfg['initial_boundary_factor']==0
    for c in case['boundary']['constraints']:assert c['value_m']==0 and not np.asarray(c['gradient']).any()
    xyz=np.asarray(mesh['nodes']);ndof=xyz.size;active=np.ones(len(mesh['tetrahedra']),bool)
    if ndof>6500:raise ValueError('independent direct reference is limited to 6500 soil DOFs')
    materials={m['id']:m for m in ground['materials']};regions={r['id']:r for r in mesh['regions']};seeds=cfg['seed']['states']
    k=np.zeros((ndof,ndof));seed_force=np.zeros(ndof);operators=[]
    for ei in range(len(active)):
        m=materials[mesh['material_ids'][ei]]['mohr_coulomb'];d=ref.elastic_tensor(m['young_kpa'],m['poisson'])
        seed=seeds[regions[mesh['region_ids'][ei]]['body']];s0=np.asarray(seed['effective_stress_compression_positive_kpa'])
        q=ref.volume(mesh,ei,independent=True);dofs=(q['nodes'][:,None]*3+np.arange(3)).reshape(-1)
        k[np.ix_(dofs,dofs)]+=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        total=np.broadcast_to(s0,(len(q['weight']),6)).copy();h=ground['initialization']['water_level_m']
        if h is not None:total[:,:3]+=np.maximum(h-q['position'][:,2],0)[:,None]*ground['initialization']['water_unit_weight_kn_m3']
        seed_force[dofs]+=-np.einsum('qai,qa,q->i',q['b'],total,q['weight'],optimize=True)
        operators.append((ref.volume(mesh,ei),d,s0,seed['void_ratio']))
    basis,error=ref.staged_constraint_space(mesh,case,{'active_elements':active,'boundary_factor':0,'displacement_m':np.zeros_like(xyz)})
    assert error==0;count=sum(b.shape[1] for b in basis.values());t=np.zeros((ndof,count));start=0
    for n,b in basis.items():t[3*n:3*n+3,start:start+b.shape[1]]=b;start+=b.shape[1]
    beams=System(mesh,case['beam_system']);state=beams.initial_state();u=np.zeros_like(xyz);factor=0.;reports=[]
    for fi,frame in enumerate(native['frames']):
        if fi:
            stage=case['stages'][fi-1]
            if stage['kind']=='boundary':factor=stage['target_factor']
            else:state=beams.change(state,u,stage)
        kb,rest,fb,eff,extf,extm,info=beams.assemble(state,factor)
        offset,binding,qb,normals=beams.space(state,u,soil_basis=t)
        # Retain beam translations, rotations and released end rotations in
        # the same direct system; no beam Schur complement is used here.
        kk=np.block([[t.T@k@t+binding.T@kb@binding,binding.T@kb@qb],
                     [qb.T@kb@binding,qb.T@kb@qb]])
        fsoil,_=physical(mesh,ground,active,ground['initialization']['surface_surcharge_kpa'],case.get('pressure_on_original_surface_only',False))
        fsoil+=factor*ref.conditions_load(mesh,case['boundary'],active)
        beam_rhs=fb-rest-kb@offset
        rhs=np.r_[t.T@(fsoil.reshape(-1)-seed_force)+binding.T@beam_rhs,qb.T@beam_rhs]
        solution=np.linalg.solve(kk,rhs);u=(t@solution[:count]).reshape(-1,3)
        beam_d=offset+binding@solution[:count]+qb@solution[count:]
        br=beams.results(state,beam_d,kb,rest,fb,eff,extf,extm,info,normals);state=br['state']
        residual=float(np.linalg.norm(kk@solution-rhs)/max(np.linalg.norm(rhs),1.));assert residual<1e-10
        actual=frame['state'];assert actual['beams']['active_members']==state['active_members']
        du=float(np.linalg.norm(np.asarray(actual['displacement_m'])-u)/max(np.linalg.norm(u),1e-12));assert du<1e-7,du
        dr=float(np.max(np.abs(np.asarray(actual['beams']['node_rotation_rad'])-state['node_rotation_rad'])))
        db=float(np.max(np.abs(np.asarray(actual['beams']['node_displacement_m'])-state['node_displacement_m'])))
        refs=float(np.max(np.abs(np.asarray(actual['beams']['reference_global_kinematics'])-state['reference_global_kinematics'])))
        assert dr<1e-9 and db<1e-9 and refs<1e-9,(dr,db,refs)
        stress=[];void=[]
        for qp,d,s0,e0 in operators:
            strain=np.einsum('qai,i->qa',qp['b'],u[qp['nodes']].reshape(-1))
            stress.extend(s0-strain@d.T);void.extend((1+e0)*np.exp(strain[:,:3].sum(axis=1))-1)
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
        reports.append({'frame':frame['name'],'soil_free_dofs':count,'retained_frame_unknowns':qb.shape[1],
                        'relative_joint_system_residual':residual,'relative_3d_displacement_error':du,
                        'maximum_beam_displacement_error_m':db,'maximum_rotation_error_rad':dr,
                        'maximum_installation_reference_error':refs,'maximum_end_action_error':error,
                        'relative_all_point_stress_error':ds,'maximum_void_ratio_error':de,
                        'global_moment_balance_kn_m':moment.tolist()})
    report={'schema':'soil-fem-independent-beam-coupling/1','reference':'full soil + frame direct NumPy solve, input seed and independently propagated installation references; all rotations and releases retained; section flexibility by virtual work; GL125 soil integration; SVD constraints',
            'frames':reports,'solver_sha256':native['provenance']['solver_binary_sha256'],
            'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (mp,gp,rp,Path(__file__),Path(ref.__file__),Path(__file__).with_name('beam_reference.py'))}}
    with output.open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();verify(args.run.resolve(),args.output.resolve())
