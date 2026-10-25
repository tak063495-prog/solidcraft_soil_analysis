"""Independent wet-interface input/flow/load/contact/force/output audit."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import json
from pathlib import Path
import meshio
import numpy as np
import fem_reference as fem
import darcy_reference
from joint_water_reference import Water
from verify_high_order import prepare,internal,physical,constraints,reduced_residual
from verify_acceptance_packet import sha

def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def audit(mesh,ground,run,exports=None):
    assert run['complete'] and run['schema']=='soil-fem-stage-result/13'
    case=run['case'];assert not case.get('beam_system') and not case.get('anchor_system'),'wet combined-structure independent assembly is a separate gate'
    assert all(s['kind'] in ('pressure','boundary','hydrostatic_water','steady_seepage','excavate') for s in case['stages']),'this wet verification gate requires explicit extension for other stage kinds'
    water=Water(mesh,ground,case,run['joint_integration_points']);q,offset=prepare(mesh,run['integration_points']);xyz=np.asarray(mesh['nodes']);positions=np.concatenate([g['position'] for g in q]);points_soil=np.asarray([water.soil[p['element_index']] for p in run['integration_points']]);regions={r['id']:r['body'] for r in mesh['regions']}
    active=np.asarray([regions[r] not in case.get('initially_inactive_bodies',[]) for r in mesh['region_ids']]);level=ground['initialization']['water_level_m'];has_water_stage=any(s['kind'] in ('hydrostatic_water','steady_seepage') for s in case['stages'])
    if level is None and has_water_stage:level=float(xyz[:,2].min())
    pressure=ground['initialization']['surface_surcharge_kpa'];factor=ground['initialization']['geostatic'].get('initial_boundary_factor',0.);had_flow=False
    pore=np.zeros(len(positions)) if level is None else np.maximum(level-positions[:,2],0)*water.gamma;pore[~points_soil]=0.;interface=water.hydrostatic(active,level)
    maximum_force=maximum_pressure=maximum_head=0.;rows=[];flow_rows=[]
    for index,frame in enumerate(run['frames']):
        state=frame['state'];previous_pore=pore.copy()
        if index:
            stage=case['stages'][index-1];kind=stage['kind']
            if kind=='excavate':
                selected=np.asarray([regions[r] in stage['bodies'] for r in mesh['region_ids']]);
                for name in stage['bodies']:
                    if name in mesh.get('excavation_sets',{}):selected[mesh['excavation_sets'][name]]=True
                active=active&~selected;pressure=stage['remaining_surface_pressure_kpa']
            elif kind=='pressure':pressure=stage['target_pressure_kpa']
            elif kind=='boundary':factor=stage['target_factor']
            if kind=='steady_seepage':
                level=stage['exterior_water_level_m'];solved=water.steady(stage['flow'],active,state.get('element_material_ids') or mesh['material_ids'],level);field=solved['field'];audit_row=next(a for a in run['flow_audits'] if a['stage_index']==index-1);native=audit_row['joint_flow'];assert native is not None
                mapping=np.asarray(native['original_node_to_hydraulic_node']);expected=solved['mapping'];assert np.array_equal(mapping[:,None]==mapping[None,:],expected[:,None]==expected[None,:])
                used=np.unique(fem.cells(mesh)[active&water.soil]);difference=float(np.max(np.abs(np.asarray(native['original_field']['head_m'])[used]-field['head'][used])));maximum_head=max(maximum_head,difference);assert difference<1e-8
                flow_rows.append(darcy_reference.verify(native['hydraulic_mesh'],native['hydraulic_case'],audit_row['solution'],state.get('element_material_ids') or mesh['material_ids']))
                for e,g in enumerate(q):
                    if not active[e]:continue
                    if not water.soil[e]:pore[offset[e]:offset[e+1]]=0.;continue
                    nodes=g['nodes'];h=field['relative_head'][nodes];psi=field['datum'][nodes[0]]-g['position'][:,2]+h[0]+g['shape']@(h-h[0]);assert psi.min()>-1e-8;pore[offset[e]:offset[e+1]]=np.maximum(psi,0)*water.gamma
                interface=water.from_heads(active,field,level);had_flow=True
            elif kind=='hydrostatic_water':
                level=stage['target_level_m'];pore=np.maximum(level-positions[:,2],0)*water.gamma;pore[~points_soil]=0.
                if had_flow:
                    for e,live in enumerate(active):
                        if not live:pore[offset[e]:offset[e+1]]=previous_pore[offset[e]:offset[e+1]]
                interface=water.hydrostatic(active,level)
            else:interface=water.retain(interface,active)
        assert np.array_equal(active,state['active_elements']);assert state['surface_pressure_kpa']==pressure and state['boundary_factor']==factor
        if has_water_stage:assert state['hydrostatic_water_level_m']==level
        assert (state.get('pore_pressure_kpa') is not None)==had_flow
        if had_flow:assert np.allclose(state['pore_pressure_kpa'],pore,rtol=1e-10,atol=1e-8)
        frame_ground={**ground,'initialization':{**ground['initialization'],'water_level_m':level}}
        force,yield_value=internal(mesh,frame_ground,state['material_states'],active,q,offset,state.get('element_material_ids'),pore)
        load,geometry=physical(mesh,frame_ground,active,pressure,case['pressure_on_original_surface_only'],state.get('element_material_ids'),water.contact)
        if case.get('boundary') and factor:load+=factor*fem.conditions_load(mesh,case['boundary'],active)
        load+=interface['load'];force+=water.audit(state,frame['joint_results'],interface)
        error=float(np.max(np.abs(load.ravel()-state['applied_nodal_load_kn'])));maximum_force=max(maximum_force,error);assert error<1e-8*max(np.linalg.norm(load),1.)
        if frame.get('joint_support'):assert water.contact.rigid_support(case,state,frame['joint_support'])['unrestrained_modes']==0
        used,fixed=constraints(mesh,active,case['settings']['kinematics']);residual=load-force
        if case.get('boundary'):
            basis,constraint_error=fem.staged_constraint_space(mesh,case,state);reduced=np.concatenate([b.T@residual[n] for n,b in basis.items()]);assert constraint_error<1e-10
            reaction=np.zeros_like(load)
            for n,b in basis.items():r=force[n]-load[n];reaction[n]=r-b@(b.T@r)
        else:
            reduced=reduced_residual(xyz,residual,fixed,used,case['settings']['kinematics']);reaction=np.where(fixed,force-load,0);reaction[np.setdiff1d(np.arange(len(xyz)),used)]=0.
        relative=float(np.linalg.norm(reduced)/max(np.linalg.norm(load),1.));assert relative<1.05*case['settings']['relative_force_tolerance']+1e-10,(frame['name'],relative)
        assert np.allclose(reaction.sum(axis=0),frame['support_reaction_kn'],rtol=1e-9,atol=1e-7);assert np.allclose(load.sum(axis=0),frame['applied_force_kn'],rtol=1e-10,atol=1e-8)
        moment=np.cross(xyz,load+reaction).sum(axis=0);assert np.allclose(moment,frame['global_moment_balance_error_kn_m'],rtol=1e-8,atol=1e-6)
        assert np.linalg.norm(state['excavation_release_load_kn'])<1e-9
        if exports:
            raw=meshio.read(exports/f'integration-points-{index:03}.vtu');ids=raw.point_data['integration_point_id'];difference=float(np.max(np.abs(raw.point_data['pore_pressure_kpa']-pore[ids])));maximum_pressure=max(maximum_pressure,difference);assert difference<1e-8
            assert np.array_equal(raw.point_data['material_stress_reference'],(~points_soil[ids]).astype(np.uint8));stress=np.asarray([state['material_states'][i]['effective_stress_compression_positive_kpa'] for i in ids]);stress[:,:3]+=pore[ids,None];assert np.allclose(raw.point_data['total_stress_compression_positive_kpa'],stress,rtol=1e-10,atol=1e-8)
        rows.append(dict(name=frame['name'],relative_force_residual=relative,physical_load_component_error_kn=error,maximum_normalized_yield=yield_value,**geometry))
    assert all(a['committed_state_unchanged_on_rejection'] for a in run['attempts'] if not a['accepted'])
    return dict(schema='soil-fem-wet-joint-independent-verification/1',passed=True,frames=rows,independent_darcy_solutions=flow_rows,maximum_physical_load_component_error_kn=maximum_force,maximum_exported_pressure_difference_kpa=maximum_pressure,maximum_independent_hydraulic_head_difference_m=maximum_head,scope='Independent input-based head constraints and dense Darcy FEM, physical pressure trace and wet-face integration, stored material/contact stress and energy/force/rigid-support reassembly, output readback. No independent general plastic mechanics or unsaturated joint PDE is claimed.')

def verify_package(mesh_path,ground_path,result_path,exports,solver=None):
    mesh=read(mesh_path);ground=fem.normalize_ground(read(ground_path));run=read(result_path);assert run['ground']==ground;case_path=result_path.parent/'stages.json';assert fem.normalize_stage_case(read(case_path))==run['case']
    result=audit(mesh,ground,run,Path(exports));export=read(Path(exports)/'export-verification.json');assert export['vtu_roundtrip']=='passed' and export['joint_point_roundtrip']=='passed' and export['result_sha256']==sha(result_path) and export['mesh_sha256']==sha(mesh_path)
    root=Path(__file__).resolve().parents[1];assert run['provenance']['solver_binary_sha256']==sha(solver or root/'target/release/soil-fem-3d.exe')
    with (result_path.parent/'stage-verification.json').open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2)
    with (result_path.parent/'accepted-stage-run.json').open('x',encoding='utf-8') as stream:json.dump(dict(schema='soil-fem-accepted-staged-run/1',solver_binary_sha256=run['provenance']['solver_binary_sha256'],inputs_and_results_sha256={str(p.resolve()):sha(p) for p in (mesh_path,ground_path,case_path,result_path)}),stream,indent=2)
    print('Wet JOINT independent audit:',len(result['frames']),'frames');return result
