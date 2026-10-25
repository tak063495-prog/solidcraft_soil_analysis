"""Independent assembly checks for general initial-stress equilibration.

This module audits stress/holding/physical loads and the retained reference.
The separate slope benchmark solves an independent linear system as well.
"""
import numpy as np
import fem_reference as ref


def verify_monitors(mesh,case,frame,previous,reference):
    xyz=np.asarray(mesh['nodes']);state=frame['state'];used=np.unique(ref.cells(mesh)[np.asarray(state['active_elements'])])
    groups={name:np.flatnonzero(np.isclose(xyz[:,2],z,atol=1e-8,rtol=0)) for name,z in case['monitors'].items()}
    groups.update({name:np.asarray(mesh['node_sets'][key]) for name,key in case.get('monitor_sets',{}).items()})
    u=np.asarray(state['displacement_m']);old=np.asarray(previous['state']['displacement_m'])
    for name,nodes in groups.items():
        nodes=np.intersect1d(nodes,used);m=frame['monitors'][name];assert m['active_node_count']==len(nodes)
        if not len(nodes):assert m['mean_displacement_z_m'] is None;continue
        z=(u-reference)[nodes,2]
        assert np.isclose(m['mean_displacement_z_m'],z.mean(),atol=1e-12,rtol=1e-10)
        assert np.isclose(m['minimum_displacement_z_m'],z.min(),atol=1e-12,rtol=1e-10)
        assert np.isclose(m['maximum_displacement_z_m'],z.max(),atol=1e-12,rtol=1e-10)
        assert np.isclose(m['stage_increment_z_m'],(u-old)[nodes,2].mean(),atol=1e-12,rtol=1e-10)


def verify_initial(mesh,ground,run,q,offset):
    from verify_high_order import internal,physical,reduced_residual,constraints
    cfg=ground['initialization']['geostatic'];g=run['geostatic'];assert g['complete'] and g['progress']==1.
    seed=g['seed'];last=g['last_committed'];all_active=np.ones(len(q),dtype=bool);case=run['case'];xyz=np.asarray(mesh['nodes'])
    inactive=set(case.get('initially_inactive_bodies',[]));regions={r['id']:r['body'] for r in mesh['regions']}
    all_active=np.asarray([regions[r] not in inactive for r in mesh['region_ids']])
    assert np.array_equal(seed['active_elements'],all_active) and not np.asarray(seed['displacement_m']).any()
    definition=cfg['seed']
    if definition['kind']=='body_states':
        regions={r['id']:r['body'] for r in mesh['regions']}
        expected=[definition['states'][regions[mesh['region_ids'][p['element_index']]]] for p in run['integration_points']]
    else:
        assert definition['mesh_sha256']==run['provenance']['mesh_sha256'];expected=[]
        for p,qinfo in zip(definition['points'],run['integration_points']):
            assert p['element_index']==qinfo['element_index'] and p['local_index']==qinfo['local_index']
            assert p['material_id']==mesh['material_ids'][p['element_index']]
            assert np.allclose(p['position_m'],qinfo['position_m'],atol=1e-8,rtol=0);expected.append(p['state'])
    assert expected==seed['material_states']
    force,_=internal(mesh,ground,seed['material_states'],all_active,q,offset)
    joint_system=None
    if case.get('joint_system'):
        from joint_reference import System as JointSystem,coupled_forces as joint_forces
        joint_system=JointSystem(mesh,case['joint_system'],run['joint_integration_points'])
    base,_=physical(mesh,ground,all_active,ground['initialization']['surface_surcharge_kpa'],case.get('pressure_on_original_surface_only',False),joint_system=joint_system)
    beam_system=None
    if case.get('beam_system'):
        from beam_reference import System,coupled_forces
        beam_system=System(mesh,case['beam_system'])
        force,seed_base,_=coupled_forces(beam_system,seed,force,base)
    else:seed_base=base
    anchor_system=None
    if case.get('anchor_system'):
        from anchor_reference import System as AnchorSystem,coupled_forces as anchor_forces
        anchor_system=AnchorSystem(mesh,case['anchor_system'])
        force,seed_base=anchor_forces(anchor_system,seed,force,seed_base)
    if joint_system is not None:force,seed_base=joint_forces(joint_system,seed,force,seed_base)
    assert np.allclose(seed['applied_nodal_load_kn'],force.reshape(-1),atol=1e-8,rtol=1e-10)
    holding=force-seed_base;assert np.allclose(g['initial_holding_load_kn'],holding.reshape(-1),atol=1e-8,rtol=1e-10)
    assert np.allclose(seed['geostatic_release_load_kn'],holding.reshape(-1),atol=1e-8,rtol=1e-10)
    assert not np.asarray(last['geostatic_release_load_kn']).any() and not np.asarray(last['excavation_release_load_kn']).any()
    assert last['boundary_factor']==cfg['initial_boundary_factor']
    physical_load=base
    if case.get('boundary') is not None:physical_load=base+cfg['initial_boundary_factor']*ref.conditions_load(mesh,case['boundary'],all_active)
    actual,_=internal(mesh,ground,last['material_states'],all_active,q,offset)
    if beam_system is not None:actual,physical_load,_=coupled_forces(beam_system,last,actual,physical_load)
    if anchor_system is not None:actual,physical_load=anchor_forces(anchor_system,last,actual,physical_load)
    if joint_system is not None:actual,physical_load=joint_forces(joint_system,last,actual,physical_load)
    assert np.allclose(last['applied_nodal_load_kn'],physical_load.reshape(-1),atol=1e-8,rtol=1e-10)
    residual=physical_load-actual
    if case.get('boundary') is not None:
        basis,error=ref.staged_constraint_space(mesh,case,last);assert error<1e-10
        reduced=np.concatenate([b.T@residual[n] for n,b in basis.items()])
    else:
        used,fixed=constraints(mesh,all_active,case['settings']['kinematics']);reduced=reduced_residual(xyz,residual,fixed,used,case['settings']['kinematics'])
    relative=np.linalg.norm(reduced)/max(np.linalg.norm(physical_load),1.)
    assert relative<min(case['settings']['relative_force_tolerance'],ground['initialization']['relative_equilibrium_tolerance'])*1.01
    assert np.isclose(relative,g['final_relative_physical_force_residual'],atol=1e-12,rtol=1e-4)
    reference=np.asarray(run['reference_displacement_m']);assert np.array_equal(reference,last['displacement_m'])
    frame=run['frames'][0];assert frame['name']=='initial_equilibrated' and frame['state']['material_states']==last['material_states']
    assert np.array_equal(frame['state']['displacement_m'],reference)
    assert frame['state']['stage_index']==frame['state']['stage_progress']==0 and not frame['state'].get('geostatic_release_load_kn')
    strain=max(float(np.max(np.abs(np.einsum('qai,i->qa',v['b'],reference[v['nodes']].reshape(-1))))) for ei,v in enumerate(q) if all_active[ei])
    assert np.isclose(strain,g['maximum_absolute_engineering_strain_from_seed'],rtol=1e-10,atol=1e-12)
    assert strain<=cfg['maximum_initial_engineering_strain']
    assert all(a['committed_state_unchanged_on_rejection'] for a in g['attempts'] if not a['accepted'])
    assert g.get('failure_reason') is None
