"""Reconstruct explicit lifecycle operations from input, before constitutive solves."""
import copy
import numpy as np
import fem_reference as ref


def verify_transition(mesh,ground,run,index,q,offset):
    from verify_high_order import internal,physical,constraints,reduced_residual
    case=run['case'];stage=case['stages'][index-1];origin=run['frames'][index-1]['state'];last=run['frames'][index]['state']
    expected=copy.deepcopy(origin);expected.setdefault('element_material_ids',mesh['material_ids'])
    structural=stage['kind'] in ('install_beams','remove_beams','install_anchors','remove_anchors')
    if structural:
        if stage['kind'] in ('install_anchors','remove_anchors'):
            from anchor_reference import System
            expected['anchors']=System(mesh,case['anchor_system']).change(origin['anchors'],np.asarray(origin['displacement_m']),stage)
            # Locate independently; coordinate inversion can differ at roundoff.
            native=next(a for a in run['transition_audits'] if a['stage']==stage['name'])['prepared_state']['anchors']
            for a,b in zip(expected['anchors']['anchors'],native['anchors']):
                assert np.allclose(a['reference_extensions_m'],b['reference_extensions_m'],atol=1e-10,rtol=1e-8)
                a['reference_extensions_m']=b['reference_extensions_m']
        else:
            from beam_reference import System
            expected['beams']=System(mesh,case['beam_system']).change(origin['beams'],np.asarray(origin['displacement_m']),stage)
        selected=[];definitions=None
    else:
        selected=[i for i,r in enumerate(mesh['region_ids']) if next(b['body'] for b in mesh['regions'] if b['id']==r) in stage['bodies']]
        assert selected and set(stage['bodies'])=={next(b['body'] for b in mesh['regions'] if b['id']==mesh['region_ids'][i]) for i in selected}
        definitions=stage['states'] if stage['kind'] in ('birth','reset_state') else stage['history'].get('states')
    if definitions is not None:assert set(definitions)==set(stage['bodies'])
    if stage['kind']=='birth' and stage.get('placement')=='cad_new_nodes':
        cells=ref.cells(mesh);existing=np.unique(cells[np.asarray(origin['active_elements'])])
        new=np.setdiff1d(np.unique(cells[selected]),existing)
        for n in new:expected['displacement_m'][n]=[0.,0.,0.]
    for ei in selected:
        body=next(b['body'] for b in mesh['regions'] if b['id']==mesh['region_ids'][ei])
        if stage['kind']=='birth':
            assert not origin['active_elements'][ei];expected['active_elements'][ei]=True
            s=definitions[body];assert not np.asarray(s['effective_stress_compression_positive_kpa']).any()
            assert not np.asarray(s['plastic_strain_compression_positive']).any()
            assert s['equivalent_plastic_strain']==s['plastic_work_kpa']==0 and s['last_return']=='elastic'
            h=ground['initialization']['water_level_m'];assert h is None or np.all(q[ei]['position'][:,2]>=h)
        elif stage['kind']=='material_set':expected['element_material_ids'][ei]=stage['material_id']
        if definitions is not None:
            expected['material_states'][offset[ei]:offset[ei+1]]=[copy.deepcopy(definitions[body]) for _ in range(offset[ei+1]-offset[ei])]
    if stage['kind']=='birth' and expected.get('display_reference_displacement_m'):
        cells=ref.cells(mesh);existing=np.unique(cells[np.asarray(origin['active_elements'])]);new=np.setdiff1d(np.unique(cells[selected]),existing)
        for n in new:expected['display_reference_displacement_m'][n]=expected['displacement_m'][n].copy()
    assert last.get('display_reference_displacement_m')==expected.get('display_reference_displacement_m')
    assert last['active_elements']==expected['active_elements'] and last['element_material_ids']==expected['element_material_ids']
    audit=next(a for a in run['transition_audits'] if a['stage']==stage['name'])
    assert audit['selected_elements']==selected and audit['prepared_state']==expected
    active=np.asarray(expected['active_elements']);ids=expected['element_material_ids']
    f,_=physical(mesh,ground,active,last['surface_pressure_kpa'],case.get('pressure_on_original_surface_only',False),ids)
    if case.get('boundary') and last.get('boundary_factor',0):f+=last['boundary_factor']*ref.conditions_load(mesh,case['boundary'],active)
    force,_=internal(mesh,ground,expected['material_states'],active,q,offset,ids)
    if case.get('beam_system'):
        from beam_reference import System,coupled_forces
        force,f,_=coupled_forces(System(mesh,case['beam_system']),expected,force,f)
    if case.get('anchor_system'):
        from anchor_reference import System,coupled_forces
        force,f=coupled_forces(System(mesh,case['anchor_system']),expected,force,f)
    holding=force-f
    used=np.unique(ref.cells(mesh)[active]);holding[np.setdiff1d(np.arange(len(mesh['nodes'])),used)]=0
    assert np.allclose(audit['initial_holding_load_kn'],holding.reshape(-1),atol=1e-7,rtol=1e-9)
    residual=f+np.asarray(audit['initial_holding_load_kn']).reshape(-1,3)-force
    if case.get('boundary'):
        basis,_=ref.staged_constraint_space(mesh,case,expected);reduced=np.concatenate([b.T@residual[n] for n,b in basis.items()])
    else:
        used,fixed=constraints(mesh,active,case['settings']['kinematics']);reduced=reduced_residual(np.asarray(mesh['nodes']),residual,fixed,used,case['settings']['kinematics'])
    relative=np.linalg.norm(reduced)/max(np.linalg.norm(f),1.)
    assert relative<1e-10 and np.isclose(relative,audit['zero_release_relative_force_residual'],atol=1e-12)
    return expected


def verify_display_reset(run,index):
    before=run['frames'][index-1]['state'];after=run['frames'][index]['state'];expected=copy.deepcopy(before)
    expected['display_reference_displacement_m']=copy.deepcopy(before['displacement_m'])
    expected['stage_index']=index;expected['stage_progress']=1.
    assert after==expected
    assert all(m['mean_displacement_z_m'] in (None,0.) for m in run['frames'][index]['monitors'].values())
