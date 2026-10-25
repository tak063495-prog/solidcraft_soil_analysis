"""CLI, CAD domain transfer, independent force/flow audits, elastic PDE, restart."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
import fem_reference as fem
import one_way_reference as ref
from verify_high_order import internal
from verify_richards import audit as flow_audit
from verify_checkpoint import invoke as execute
from verify_one_way_cases import normalized
from verify_acceptance_packet import sha
from export_one_way_construction import export

def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def invoke(solver,args,log,success=True):execute([str(solver)]+[str(a) for a in args],log,0 if success else 1)
def write(path,data):
    with path.open('x',encoding='utf-8') as stream:json.dump(data,stream,indent=2)

def domain(mesh,active):
    out=copy.deepcopy(mesh);owners=fem.topology(mesh)
    exposed={key for key,items in owners.items() if sum(active[e] for e in items)==1}
    out['face_sets']={name:[f for f in fs if tuple(sorted(f)) in exposed] for name,fs in mesh['face_sets'].items()}
    out['face_sets']={n:fs for n,fs in out['face_sets'].items() if fs}
    return out

def audit(mesh,run):
    assert run['complete'];ground=run['ground'];case=run['case']['one_way'];geos=[ref.geometry(mesh,e) for e in range(len(mesh['tetrahedra']))]
    mechanics=[fem.volume(mesh,e) for e in range(len(geos))];offset=np.cumsum([0]+[len(g['weight']) for g in mechanics]).tolist()
    maximum_force=maximum_residual=maximum_water_error=0.;rows=[];segments=[];initial=run['initialization']['last_committed'];segment=[dict(name='initial',state=initial['hydraulic'],balance=None)]
    active_previous=initial['mechanical']['active_elements'];previous=initial
    for frame in run['frames']:
        state=frame['state'];active=state['mechanical']['active_elements'];projected=domain(mesh,active)
        if active!=active_previous:
            segments.append((active_previous,segment));segment=[dict(name='boundary_jump',state=state['hydraulic'],balance=None)]
            record=next(r for r in run['excavations'] if r['name']==frame['name'])
            assert record['origin_hydraulic']==previous['hydraulic']
            removal=record['removal'];mass=0.
            for cell in removal['removed_cells']:
                e=cell['element_index'];g=geos[e];p=next(m['retention'] for m in case['flow']['materials'] if m['id']==mesh['material_ids'][e]);h=previous['hydraulic']['relative_cell_head_m'][e]
                psi=case['flow']['head_datum_m']-g['position'][:,2]+h;theta,_,_=ref.retention(psi,p)
                water=float(g['weight']@(theta+p['saturated_specific_storage_per_m']*np.maximum(psi,0)))
                assert abs(water-cell['water_equivalent_volume_m3'])<1e-11*max(water,1.)
                assert cell['relative_head_m']==h and cell['storage_offset_m3']==previous['hydraulic']['cell_storage_offset_m3'][e]
                assert cell['outflow_m3_s']==previous['hydraulic']['cell_outflow_m3_s'][e];mass+=water
            assert abs(mass-removal['water_removed_m3'])<1e-11*max(mass,1.)
            for i,j in enumerate(removal['new_to_old_face']):
                assert removal['new_faces'][i]['corners']==removal['old_faces'][j]['corners']
                assert removal['state']['relative_face_head_m'][i]==previous['hydraulic']['relative_face_head_m'][j]
            for e,live in enumerate(active):
                if live:
                    for key in ('relative_cell_head_m','cell_storage_offset_m3','cell_outflow_m3_s'):assert state['hydraulic'][key][e]==previous['hydraulic'][key][e]
            for i in record['removed_mechanical_points']:assert previous['mechanical']['material_states'][i]==run['last_committed']['mechanical']['material_states'][i]
        if frame['phase']=='excavation_release':
            segment[0]=dict(name='boundary_jump',state=state['hydraulic'],balance=None)
            assert frame['hydraulic_balance'] is None
        else:segment.append(dict(name=frame['name'],state=state['hydraulic'],balance=frame['hydraulic_balance']))
        reservoirs=next(p['reservoirs'] for p in case['exterior_periods'] if p['name']==frame['name']) if frame['phase']=='transient' else next(p['reservoirs'] for p in case['exterior_periods'] if p['name']==next(e['before_period'] for e in run['case']['excavations'] if e['name']==frame['name']))
        base=ref.base_load(mesh,ground,case['mechanics'],active)
        water=ref.fluid_load(projected,ground,case,state['hydraulic'],reservoirs,geos,base)
        for key in ('fluid_internal_nodal_force_kn','water_body_nodal_load_kn','exterior_nodal_load_kn','physical_nodal_load_kn','equivalent_mechanical_rhs_kn'):
            delta=float(np.max(np.abs(water[key]-state['fluid'][key])));maximum_force=max(maximum_force,delta)
            assert delta<1e-8*max(np.linalg.norm(water[key]),1.),(frame['name'],key,delta)
        assert abs(water['water_equivalent_mass_volume_m3']-state['fluid']['water_equivalent_mass_volume_m3'])<1e-11*max(water['water_equivalent_mass_volume_m3'],1.)
        force,_=internal(mesh,ground,state['mechanical']['material_states'],active,mechanics,offset,pore_pressure_kpa=[0.]*len(run['mechanical_integration_points']))
        q,affine=ref.constraint_space(mesh,case['mechanics'],ground,active)
        rhs=water['equivalent_mechanical_rhs_kn']+np.asarray(frame['holding_load_kn']);relative=float(np.linalg.norm(q.T@(rhs-force.ravel()))/max(np.linalg.norm(rhs),1.));maximum_residual=max(maximum_residual,relative)
        assert relative<1.05*case['mechanics']['settings']['relative_force_tolerance']+1e-10,(frame['name'],relative)
        error=water['water_equivalent_mass_volume_m3']+frame['cumulative_removed_water_m3']-run['initial_water_m3']-frame['cumulative_net_water_inflow_m3'];maximum_water_error=max(maximum_water_error,abs(error))
        assert abs(error-frame['signed_global_water_budget_error_m3'])<1e-10*max(run['initial_water_m3'],1.)
        assert abs(error)<1e-8*max(run['initial_water_m3'],1.)
        rows.append(dict(name=frame['name'],phase=frame['phase'],time_s=state['hydraulic']['time_s'],relative_mechanical_residual=relative,water_budget_error_m3=error))
        active_previous=active;previous=state
    segments.append((active_previous,segment));water_checks=[]
    for active,frames in segments:
        if len(frames)==1:continue
        for frame in frames:frame['points']=[]
        times={f['state']['time_s'] for f in frames[1:]}
        attempts=[dict(period=a['name'],accepted=True,time_step_s=a['dt_s'],start_time_s=a['time_s']) for a in run['attempts'] if a['accepted'] and a['phase']=='transient' and a['time_s']+a['dt_s'] in times]
        representative=next(f for f in run['frames'] if f['state']['mechanical']['active_elements']==active)
        water_run=dict(schema='soil-fem-richards-result/1',complete=True,active_elements=active,cell_element_indices=np.flatnonzero(active).tolist(),material_ids=mesh['material_ids'],case=case['flow'],faces=representative['hydraulic_faces'],frames=frames,attempts=attempts)
        water_checks.append(flow_audit(domain(mesh,active),water_run))
    return dict(schema='soil-fem-construction-independent-reassembly/1',passed=True,maximum_force_component_difference_kn=maximum_force,maximum_relative_mechanical_residual=maximum_residual,maximum_global_water_budget_error_m3=maximum_water_error,frames=rows,hydraulic_segments=water_checks,scope='Independent geometry, transfer, removed mass, stored-stress force and local RT0 balance reassembly; no independent plastic 3D solution is inferred.')

def elastic(mesh,run):
    assert all(not any(s['plastic_strain_compression_positive']) for s in run['last_committed']['mechanical']['material_states'])
    case=run['case']['one_way'];ends=[run['frames'][i] for i in range(len(run['frames'])) if run['frames'][i]['phase']=='transient' and (i+1==len(run['frames']) or run['frames'][i+1]['phase']=='excavation_release')];rows=[]
    for frame in ends:
        active=frame['state']['mechanical']['active_elements'];projected=domain(mesh,active)
        reservoirs=next(p['reservoirs'] for p in case['exterior_periods'] if p['name']==frame['name'])
        solved=[ref.ElasticReference(projected,run['ground'],case,order,active).solve(frame['state']['hydraulic'],reservoirs)['displacement_m'] for order in (5,7,9)]
        nodes=np.unique(fem.cells(mesh)[np.asarray(active,dtype=bool)]);native=np.asarray(frame['state']['mechanical']['displacement_m'])[nodes]
        relative=float(np.linalg.norm(native-solved[2][nodes])/np.linalg.norm(solved[2][nodes]));quadrature=float(np.linalg.norm(solved[0][nodes]-solved[1][nodes])/np.linalg.norm(solved[1][nodes]));refined=float(np.linalg.norm(solved[1][nodes]-solved[2][nodes])/np.linalg.norm(solved[2][nodes]))
        assert relative<1e-4 and quadrature<1e-5 and refined<max(quadrature*.6,1e-11) and refined<1e-7,(frame['name'],relative,quadrature,refined)
        rows.append(dict(name=frame['name'],relative_displacement_l2_difference=relative,independent_quadrature_5_to_7_relative_difference=quadrature,independent_quadrature_7_to_9_relative_difference=refined,maximum_displacement_component_difference_m=float(np.max(np.abs(native-solved[2][nodes])))))
    return dict(passed=True,comparisons=rows,scope='Independent 3D elastic PDE on each actual active domain, given hydraulic heads. No native displacement, stress update or tangent generates the reference.')

def profile(solver,mesh,ground,case,folder,reference=False):
    folder.mkdir();write(folder/'mesh.json',mesh);write(folder/'ground.json',ground);write(folder/'case.json',case)
    invoke(solver,['solve-one-way-construction',folder/'mesh.json',folder/'ground.json',folder/'case.json',folder/'result.json'],folder/'solve.log')
    run=read(folder/'result.json');checked=audit(mesh,run);write(folder/'verification.json',checked)
    if reference:write(folder/'elastic-comparison.json',elastic(mesh,run))
    export(folder/'mesh.json',folder/'result.json',solver,folder/'results')
    return run,checked

def restart(solver,folder,run,counts):
    rows=[]
    for count in counts:
        cp=folder/f'checkpoint-{count}.json';out=folder/f'resumed-{count}.json'
        invoke(solver,['checkpoint-one-way-construction',folder/'mesh.json',folder/'ground.json',folder/'case.json',count,cp],folder/f'checkpoint-{count}.log')
        invoke(solver,['resume-one-way-construction',folder/'mesh.json',folder/'ground.json',folder/'case.json',cp,out],folder/f'resume-{count}.log')
        assert normalized(read(out))==normalized(run)
        rows.append(dict(profile=folder.name,count=count,all_states_frames_attempts_and_archives_exact=True))
    return rows

def verify(solver,output,cad=None):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False);checks={};restarts=[];elastic_checks={}
    for model in ('mc','mcc'):
        for order in ('tet4','tet10'):
            label=f'{model}-{order}';folder=output/label
            mesh=read(root/f'examples/seepage-column-{order}-mesh.json');ground=read(root/f'examples/ground-one-way-{model}.json');case=read(root/f'examples/one-way-construction-{model}.json')
            run,checked=profile(solver,mesh,ground,case,folder,model=='mc');checks[label]=checked
            restarts+=restart(solver,folder,run,[1,2,3,4,5])
            if model=='mc':elastic_checks[label]=read(folder/'elastic-comparison.json')
    mesh=read(root/'examples/seepage-column-tet4-mesh.json');ground=read(root/'examples/ground-one-way-mcc.json');case=read(root/'examples/one-way-construction-mcc.json')
    case['verification_reject_once_excavations']=['cut'];run,checked=profile(solver,mesh,ground,case,output/'cutback');checks['cutback']=checked
    assert any(not a['accepted'] and a['candidate_changed'] for a in run['attempts'])
    assert all(a['committed_pair_sha256_before']==a['committed_pair_sha256_after'] for a in run['attempts'] if not a['accepted'])
    clean=copy.deepcopy(case);clean['verification_reject_once_excavations']=[];clean['excavations'][0]['initial_increment']=0.25
    direct,_=profile(solver,mesh,ground,clean,output/'clean-release')
    assert [f['state'] for f in run['frames']]==[f['state'] for f in direct['frames']]
    restarts+=restart(solver,output/'cutback',run,[3])
    case=read(root/'examples/one-way-construction-mc.json');case['one_way']['store_transfer_integration_points']=False
    omitted,checked=profile(solver,mesh,read(root/'examples/ground-one-way-mc.json'),case,output/'omitted');checks['omitted']=checked
    ordinary=read(output/'mc-tet4/result.json')
    def primary(s):
        s=copy.deepcopy(s);s['fluid']['hydraulic_points']=[];s['fluid']['mechanical_points']=[];return s
    assert [primary(f['state']) for f in ordinary['frames']]==[primary(f['state']) for f in omitted['frames']]
    restarts+=restart(solver,output/'omitted',omitted,[3])
    if cad:
        for model in ('mc','mcc'):
            mesh=read(cad/'mesh.json');ground=read(root/'examples/ground-local.json');case=read(root/f'examples/one-way-construction-{model}.json')
            for material in ground['materials']:
                material['unit_weight_kn_m3']=16.
                if model=='mc':
                    material.pop('mcc');material['mohr_coulomb']=dict(young_kpa=10000.,poisson=0.25,cohesion_kpa=500.,friction_deg=30.,dilation_deg=0.,tension_cutoff_kpa=100.,flow='non_associated')
            if model=='mc':
                for layer in ground['layers']:layer['consolidation']=dict(definition='not_applicable');layer['void_ratio']=dict(mode='explicit',value=0.7)
            base=case['one_way'];base['flow']['initial_head']['head_m']=0.;base['initial_porosity_tolerance']=0.08
            base['mechanics']['settings']['kinematics']='side_rollers';base['mechanics']['pressure_on_original_surface_only']=True
            base['mechanics']['monitors']={};base['mechanics']['monitor_sets']={'floor':'Pit:floor','surface':'Pit:surrounding_surface'}
            for material in base['flow']['materials']:
                material['saturated_conductivity_m_s']=(np.eye(3)*1e-5).tolist();material['retention']['alpha_per_m']=0.15
                if model=='mcc':material['retention']['saturated_water_content']=0.55
            for i,p in enumerate(base['flow']['analysis']['periods']):
                p['end_time_s']=100.*(i+1);p['boundaries'][1]['face_set']='free_surface'
                if i:p['boundaries'] += [dict(face_set=f'Pit:{name}',condition=dict(kind='atmospheric_no_ponding',rainfall_vertical_m_s=1e-6)) for name in ('floor','wall')]
            case['excavations'][0]['bodies']=['Pit'];folder=output/f'cad-{model}'
            run,checked=profile(solver,mesh,ground,case,folder,model=='mc');checks[f'cad-{model}']=checked
            restarts+=restart(solver,folder,run,[1,2,3])
            write(folder/'cad-input-link.json',dict(mesh_contents_equal=mesh==read(cad/'mesh.json'),original_mesh_sha256=sha(cad/'mesh.json'),analysis_mesh_sha256=sha(folder/'mesh.json'),source_ground_sha256=sha(root/'examples/ground-local.json'),analysis_ground_sha256=sha(folder/'ground.json'),property_derivation='Dry bulk weight 16 kN/m3, explicitly assigned water retention; MC uses a strong elastic verification material. Original CAD cells and Pit Boolean membership unchanged.'))
            if model=='mc':elastic_checks['cad-mc']=read(folder/'elastic-comparison.json')
    # Invalid future boundary must produce no accepted result file.
    bad=read(output/'mc-tet4/case.json');bad['one_way']['flow']['analysis']['periods'][1]['boundaries'].pop();write(output/'invalid-case.json',bad)
    invoke(solver,['solve-one-way-construction',output/'mc-tet4/mesh.json',output/'mc-tet4/ground.json',output/'invalid-case.json',output/'invalid-result.json'],output/'invalid.log',False)
    assert not (output/'invalid-result.json').exists()
    folder=output/'mc-tet4';before=sha(folder/'result.json')
    invoke(solver,['solve-one-way-construction',folder/'mesh.json',folder/'ground.json',folder/'case.json',folder/'result.json'],output/'existing-output.log',False)
    assert sha(folder/'result.json')==before
    failing=read(folder/'case.json');failing['one_way']['flow']['settings']['maximum_iterations']=1;failing['one_way']['flow']['settings']['minimum_time_step_s']=100.
    failing['one_way']['flow']['analysis']['periods'][1]['boundaries'][1]['condition']=dict(kind='atmospheric_no_ponding',rainfall_vertical_m_s=1e-3)
    write(output/'failing-case.json',failing)
    invoke(solver,['solve-one-way-construction',folder/'mesh.json',folder/'ground.json',output/'failing-case.json',output/'failing-result.json'],output/'genuine-exhaustion.log',False)
    failure=read(output/'failing-result.json');assert not failure['complete'] and 'iteration' in failure['failure_reason']
    assert failure['last_committed']==failure['frames'][-1]['state'] and failure['last_committed']['hydraulic']['time_s']==200. and len(failure['excavations'])==1
    assert failure['attempts'][-1]['committed_pair_sha256_before']==failure['attempts'][-1]['committed_pair_sha256_after']
    report=dict(schema='soil-fem-one-way-construction-cli-verification/1',passed=True,solver_sha256=sha(solver),cad_exercised=bool(cad),restarts=restarts,profiles=checks,elastic_comparisons=elastic_checks,cutback_clean_release_exact=True,omitted_primary_exact=True,invalid_future_boundary_rejected=True,existing_output_preserved=True,genuine_post_cut_exhaustion_preserves_committed_domain_and_archive=True,maximum_independent_mechanical_residual=max(c['maximum_relative_mechanical_residual'] for c in checks.values()),maximum_water_budget_error_m3=max(c['maximum_global_water_budget_error_m3'] for c in checks.values()))
    write(output/'verification.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('profiles','restarts','elastic_comparisons')},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--cad',type=Path);a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve(),a.cad.resolve() if a.cad else None)
