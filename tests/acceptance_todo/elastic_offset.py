"""NEW bounded reference: elastic MC soil + offset plate/beam + hydrostatics.

Source contracts: todo.ja.md B-02; interim/2026-11-14-plate6-curved-offset/
README.ja.md, "検証契約と証拠". Independent GL125 soil, section-flexibility
frame, explicit plate/beam DOFs and independently propagated installation
references. Does not extend any existing verifier's scope or consume saved
solver stress/strain as an input. No plastic soil/JOINT/anchor/BIRTH claims.
"""
import copy
import numpy as np
import fem_reference as fem
from plate6_offset_reference import System
from verify_high_order import physical


def elastic_offset_gates(mesh, ground, case, run):
    ground=fem.normalize_ground(copy.deepcopy(ground))
    case=fem.normalize_stage_case(copy.deepcopy(case))
    assert not case.get('joint_system') and not case.get('anchor_system')
    assert not case.get('initially_inactive_bodies') and run['complete']
    assert all(s['kind'] in {'install_beams','remove_beams','boundary','hydrostatic_water'} for s in case['stages'])
    seed=ground['initialization']['geostatic']['seed']['states']
    assert ground['initialization']['geostatic']['initial_boundary_factor']==0.
    assert all(c['value_m']==0. and not np.asarray(c.get('gradient',[0.,0.,0.])).any() for c in case['boundary']['constraints'])
    xyz=np.asarray(mesh['nodes']); nd=xyz.size
    assert nd<=2000, 'small independent dense reference only'
    active=np.ones(len(fem.cells(mesh)),bool)
    materials={m['id']:m['mohr_coulomb'] for m in ground['materials']}
    bodies={r['id']:r['body'] for r in mesh['regions']}
    k=np.zeros((nd,nd));operators=[]
    for e in range(len(active)):
        material=materials[mesh['material_ids'][e]]
        d=fem.elastic_tensor(material['young_kpa'],material['poisson'])
        initial=seed[bodies[mesh['region_ids'][e]]]
        assert initial['equivalent_plastic_strain']==initial['plastic_work_kpa']==0.
        q=fem.volume(mesh,e,independent=True,independent_order=5)
        ids=(3*q['nodes'][:,None]+np.arange(3)).ravel()
        k[np.ix_(ids,ids)]+=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        operators.append((q,fem.volume(mesh,e),d,ids,initial))
    basis,error=fem.staged_constraint_space(mesh,case,dict(active_elements=active,boundary_factor=0.,displacement_m=np.zeros_like(xyz)))
    assert error==0.
    count=sum(b.shape[1] for b in basis.values());t=np.zeros((nd,count));start=0
    for n,b in basis.items():t[3*n:3*n+3,start:start+b.shape[1]]=b;start+=b.shape[1]
    model=System(mesh,case['beam_system']);history=model.initial_state()
    u=np.zeros_like(xyz);factor=0.;level=ground['initialization']['water_level_m']
    assert len(run['frames'])==len(case['stages'])+1
    for i,frame in enumerate(run['frames']):
        if i:
            stage=case['stages'][i-1]
            if stage['kind']=='boundary':factor=stage['target_factor']
            elif stage['kind']=='hydrostatic_water':level=stage['target_level_m']
            else:history=model.change(history,u,stage)
        current=copy.deepcopy(ground);current['initialization']['water_level_m']=level
        seed_force=np.zeros(nd)
        for q,_,_,ids,initial in operators:
            total=np.tile(initial['effective_stress_compression_positive_kpa'],(len(q['weight']),1))
            if level is not None:
                total[:,:3]+=ground['initialization']['water_unit_weight_kn_m3']*np.maximum(level-q['position'][:,2],0.)[:,None]
            seed_force[ids]-=np.einsum('qai,qa,q->i',q['b'],total,q['weight'])
        kb,rest,fb,eff,extf,extm,info=model.assemble(history,factor)
        offset,binding,qb,normals=model.space(history,u,soil_basis=t)
        matrix=np.block([[t.T@k@t+binding.T@kb@binding,binding.T@kb@qb],
                         [qb.T@kb@binding,qb.T@kb@qb]])
        load,_=physical(mesh,current,active,ground['initialization']['surface_surcharge_kpa'],case['pressure_on_original_surface_only'])
        load+=factor*fem.conditions_load(mesh,case['boundary'],active)
        brhs=fb-rest-kb@offset
        rhs=np.r_[t.T@(load.ravel()-seed_force)+binding.T@brhs,qb.T@brhs]
        solution=np.linalg.solve(matrix,rhs)
        assert np.linalg.norm(matrix@solution-rhs)/max(np.linalg.norm(rhs),1.)<=1e-10
        u=(t@solution[:count]).reshape(-1,3)
        br=model.results(history,offset+binding@solution[:count]+qb@solution[count:],kb,rest,fb,eff,extf,extm,info,normals)
        history=br['state'];actual=frame['state']
        assert np.linalg.norm(np.asarray(actual['displacement_m'])-u)/max(np.linalg.norm(u),1e-30)<=1e-7
        for key in ('node_displacement_m','node_rotation_rad','reference_global_kinematics','plate_reference_global_kinematics'):
            np.testing.assert_allclose(actual['beams'][key],history[key],atol=1e-9,rtol=0.)
        assert actual['beams']['active_members']==history['active_members']
        assert actual['beams']['active_plates']==history['active_plates']
        stresses=[];void=[]
        for _,q,d,_,initial in operators:
            strain=np.einsum('qai,i->qa',q['b'],u[q['nodes']].ravel())
            stresses.extend(np.asarray(initial['effective_stress_compression_positive_kpa'])-strain@d.T)
            void.extend((1+initial['void_ratio'])*np.exp(strain[:,:3].sum(axis=1))-1)
        points=actual['material_states']
        assert all(p['equivalent_plastic_strain']==p['plastic_work_kpa']==0. for p in points)
        stress=np.asarray(stresses)
        assert np.linalg.norm(np.asarray([p['effective_stress_compression_positive_kpa'] for p in points])-stress)/max(np.linalg.norm(stress),1e-30)<=1e-7
        np.testing.assert_allclose([p['void_ratio'] for p in points],void,atol=1e-8,rtol=0.)
        assert {p['id'] for p in frame['beam_results'].get('plates',[])}==set(br.get('plates',{}))
        for plate in frame['beam_results'].get('plates',[]):
            reference=br['plates'][plate['id']]
            for key in ('membrane_energy_kn_m','bending_energy_kn_m','transverse_shear_energy_kn_m','drilling_penalty_energy_kn_m'):
                np.testing.assert_allclose(plate[key],reference[key],atol=1e-8,rtol=1e-7)
            for actual_point,reference_point in zip(plate['samples'],reference['samples'],strict=True):
                for key in reference_point:
                    np.testing.assert_allclose(actual_point[key],reference_point[key],atol=1e-8,rtol=1e-7)
