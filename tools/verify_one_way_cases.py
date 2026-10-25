"""Fresh-process ONE_WAY frames, force audits, elastic references and restart."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from one_way_reference import audit,ElasticReference,fluid_reservoirs
from export_one_way import export
from verify_checkpoint import invoke
from verify_acceptance_packet import sha


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,v):
    with p.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
def normalized(v):
    if isinstance(v,dict):return {k:(0. if k=='elapsed_seconds' else normalized(x)) for k,x in v.items()}
    if isinstance(v,list):return [normalized(x) for x in v]
    return v


def fixture(solver,mesh,ground,case,folder,with_export=True):
    folder.mkdir(parents=True,exist_ok=False)
    for n,v in [('mesh',mesh),('ground',ground),('case',case)]:write(folder/f'{n}.json',v)
    invoke([str(solver),'solve-one-way',str(folder/'mesh.json'),str(folder/'ground.json'),str(folder/'case.json'),str(folder/'result.json')],folder/'solve.log')
    run=read(folder/'result.json');assert run['complete'] and run['solver_binary_sha256']==sha(solver)
    checked=audit(mesh,run);write(folder/'verification.json',checked)
    if with_export:export(folder/'mesh.json',folder/'result.json',solver,folder/'results')
    return run,checked


def elastic_check(mesh,run,folder,orders=(5,7)):
    # Hydraulic heads are inputs. No native mechanical displacement/stress is
    # supplied to either independent elastic system.
    references={order:ElasticReference(mesh,run['ground'],run['case'],order) for order in orders}
    selected=[next(f for f in run['frames'] if f['initialization_progress']==1.)]
    for p in run['case']['flow']['analysis']['periods']:
        selected.append(next(f for f in reversed(run['frames']) if f['name']==p['name'] and f['phase']=='transient'))
    rows=[]
    for i,f in enumerate(selected):
        assert all(not np.any(s['plastic_strain_compression_positive']) for s in f['state']['mechanical']['material_states'])
        values={order:r.solve(f['state']['hydraulic'],fluid_reservoirs(run,f)) for order,r in references.items()}
        finest=values[orders[-1]]['displacement_m'];native=np.asarray(f['state']['mechanical']['displacement_m']);denom=max(np.linalg.norm(finest),1e-30)
        relative=float(np.linalg.norm(native-finest)/denom);maximum=float(np.max(np.abs(native-finest)));refinement=float(np.linalg.norm(values[orders[0]]['displacement_m']-finest)/denom)
        assert relative<.005 and maximum<1e-5,(f['name'],relative,maximum)
        assert refinement<.002,(f['name'],refinement)
        rows.append(dict(name=f['name'],time_s=f['state']['hydraulic']['time_s'],displacement_relative_l2_difference=relative,maximum_displacement_component_difference_m=maximum,reference_quadrature_relative_difference=refinement,true_linear_residual=values[orders[-1]]['relative_linear_residual']))
        write(folder/f'elastic-reference-{i}.json',dict(name=f['name'],input_hydraulic_state=f['state']['hydraulic'],displacement_m=finest.tolist(),scope='Independent 3D elastic mechanical PDE given hydraulic primary DOFs; not an independent coupled hydraulic/plastic solution.'))
    report=dict(schema='soil-fem-one-way-elastic-comparison/1',passed=True,quadrature_orders=list(orders),comparisons=rows)
    write(folder/'elastic-verification.json',report);return report


def restart(solver,folder,run,counts):
    rows=[];args=[str(folder/f'{n}.json') for n in ('mesh','ground','case')]
    for i,count in enumerate(counts):
        cp=folder/f'checkpoint-{i}.json';out=folder/f'resumed-{i}.json'
        invoke([str(solver),'checkpoint-one-way',*args,str(count),str(cp)],folder/f'checkpoint-{i}.log')
        saved=read(cp);assert not saved['run']['complete'] and saved['run']['runtime']
        invoke([str(solver),'resume-one-way',*args,str(cp),str(out)],folder/f'resumed-{i}.log')
        assert normalized(read(out))==normalized(run)
        rows.append(dict(folder=folder.name,count=count,time_s=saved['run']['last_committed']['hydraulic']['time_s'],initialization_progress=saved['run']['runtime']['initialization_progress'],all_states_frames_and_attempts_exact=True))
    bad=read(folder/'checkpoint-0.json');bad['run']['last_committed']['fluid']['water_body_nodal_load_kn'][2]+=1.;write(folder/'bad-checkpoint.json',bad)
    invoke([str(solver),'resume-one-way',*args,str(folder/'bad-checkpoint.json'),str(folder/'bad-result.json')],folder/'bad-checkpoint.log',1);assert not (folder/'bad-result.json').exists()
    changed=read(folder/'case.json');changed['bishop']={'kind':'constant','value':0.};write(folder/'changed-case.json',changed)
    invoke([str(solver),'resume-one-way',args[0],args[1],str(folder/'changed-case.json'),str(folder/'checkpoint-0.json'),str(folder/'changed-result.json')],folder/'changed-case.log',1);assert not (folder/'changed-result.json').exists()
    return rows


def verify(solver,output,cad=None):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False);profiles={};checks={};references={};restarts=[]
    for model in ('mc','mcc'):
        for order in ('tet4','tet10'):
            label=f'{model}-{order}';folder=output/label;mesh=read(root/f'examples/seepage-column-{order}-mesh.json');ground=read(root/f'examples/ground-one-way-{model}.json');case=read(root/f'examples/one-way-{model}.json')
            run,checked=fixture(solver,mesh,ground,case,folder);profiles[label]=run;checks[label]=checked
            assert any(a['phase']=='initialization' and a['accepted'] for a in run['attempts'])
            if model=='mc':references[label]=elastic_check(mesh,run,folder)
            else:assert any(a['phase']=='transient' and a['plastic_candidate_changed'] for a in run['attempts'])
            periods=case['flow']['analysis']['periods'];indices=[1,2,3]+[next(i for i,f in enumerate(run['frames']) if f['phase']=='transient' and f['state']['hydraulic']['time_s']==p['end_time_s']) for p in periods[:2]]+[len(run['frames'])-2]
            restarts+=restart(solver,folder,run,indices)
    mesh=read(root/'examples/seepage-column-tet4-mesh.json');ground=read(root/'examples/ground-one-way-mc.json');case=read(root/'examples/one-way-mc.json')
    case['store_transfer_integration_points']=False;case['flow']['analysis']['periods']=case['flow']['analysis']['periods'][:1];case['flow']['analysis']['periods'][0]['end_time_s']=100.;case['exterior_periods']=case['exterior_periods'][:1]
    omitted,checked=fixture(solver,mesh,ground,case,output/'omitted-points');profiles['omitted-points']=omitted;checks['omitted-points']=checked
    stored=copy.deepcopy(next(f['state'] for f in profiles['mc-tet4']['frames'] if f['phase']=='transient' and f['state']['hydraulic']['time_s']==100.))
    stored['fluid']['hydraulic_points']=[];stored['fluid']['mechanical_points']=[];assert omitted['last_committed']==stored
    restarts+=restart(solver,output/'omitted-points',omitted,[2])
    mesh=read(root/'examples/seepage-column-tet10-mesh.json');ground=read(root/'examples/ground-one-way-mcc.json');case=read(root/'examples/one-way-mcc.json')
    case['flow']['analysis']['periods']=case['flow']['analysis']['periods'][:1];case['flow']['analysis']['periods'][0]['end_time_s']=100.;case['exterior_periods']=case['exterior_periods'][:1]
    injected=copy.deepcopy(case);injected['verification_reject_once_at_times_s']=[100.]
    cut,checked=fixture(solver,mesh,ground,injected,output/'cutback',False)
    rejected=[a for a in cut['attempts'] if not a['accepted']];assert any(a['hydraulic_candidate_changed'] and a['mechanical_candidate_changed'] and a['plastic_candidate_changed'] for a in rejected)
    assert all(a['committed_pair_sha256_before']==a['committed_pair_sha256_after'] for a in rejected)
    clean=copy.deepcopy(case);bc=case['flow']['analysis']['periods'][0]['boundaries'];clean['flow']['analysis']['periods']=[dict(name=f'clean-{i}',end_time_s=f['state']['hydraulic']['time_s'],boundaries=bc) for i,f in enumerate(cut['frames']) if f['phase']=='transient']
    clean['exterior_periods']=[dict(name=p['name'],reservoirs=[]) for p in clean['flow']['analysis']['periods']]
    direct,_=fixture(solver,mesh,ground,clean,output/'clean',False);assert [f['state'] for f in cut['frames']]==[f['state'] for f in direct['frames']]
    restarts+=restart(solver,output/'cutback',cut,[3])
    exhausted=copy.deepcopy(case);exhausted['mechanics']['settings']['maximum_newton_iterations']=4;exhausted['flow']['settings']['minimum_time_step_s']=100.
    write(output/'exhaustion-case.json',exhausted)
    args=[str(solver),'solve-one-way',str(output/'mcc-tet10/mesh.json'),str(output/'mcc-tet10/ground.json')]
    invoke([*args,str(output/'exhaustion-case.json'),str(output/'exhaustion-result.json')],output/'exhaustion.log',1)
    failed=read(output/'exhaustion-result.json');a=failed['attempts'][-1];assert not failed['complete'] and a['hydraulic_candidate_changed'] and not a['accepted'] and 'Newton iteration limit' in a['reason']
    assert a['committed_pair_sha256_before']==a['committed_pair_sha256_after'] and failed['last_committed']==failed['frames'][-1]['state'] and failed['last_committed']['hydraulic']['time_s']==0.
    invalid=copy.deepcopy(case);invalid['exterior_periods'][0]['reservoirs']=[dict(face_set='missing',level_m=0.)];write(output/'invalid-case.json',invalid)
    invoke([*args,str(output/'invalid-case.json'),str(output/'invalid-result.json')],output/'invalid.log',1);assert not (output/'invalid-result.json').exists()
    original=sha(output/'mcc-tet10/result.json');invoke([*args,str(output/'mcc-tet10/case.json'),str(output/'mcc-tet10/result.json')],output/'existing-output.log',1);assert sha(output/'mcc-tet10/result.json')==original
    if cad:
        mesh=read(cad/'mesh.json');ground=read(cad/'ground.json');case=read(root/'examples/one-way-mc.json')
        for m in ground['materials']:m['unit_weight_kn_m3']=16.;m['mohr_coulomb']['cohesion_kpa']=500.
        ground['initialization']['geostatic']['verification_reject_once']=False
        case['mechanics']=read(root/'examples/stages-geostatic-slope.json');case['mechanics']['stages']=[];case['mechanics']['verification_reject_once_stages']=[]
        law=case['flow']['materials'][0];law['retention']['alpha_per_m']=.15;law['saturated_conductivity_m_s']=(np.eye(3)*1e-4).tolist()
        assert len(ground['materials'])==1,'slope reference fixture has one material body'
        case['flow']['materials']=[copy.deepcopy(law)]
        case['flow']['materials'][0]['id']=ground['materials'][0]['id']
        case['flow']['analysis']['periods']=case['flow']['analysis']['periods'][:2]
        for i,p in enumerate(case['flow']['analysis']['periods']):
            p['end_time_s']=(i+1)*1000.;p['boundaries'][1]['face_set']='free_surface';p['boundaries'][1]['condition']['rainfall_vertical_m_s']=5e-5 if i==0 else 0.
        case['flow']['settings'].update(initial_time_step_s=1000.,maximum_time_step_s=1000.)
        case['exterior_periods']=[dict(name=p['name'],reservoirs=[]) for p in case['flow']['analysis']['periods']]
        run,checked=fixture(solver,mesh,ground,case,output/'cad-slope');profiles['cad-slope']=run;checks['cad-slope']=checked;references['cad-slope']=elastic_check(mesh,run,output/'cad-slope')
        assert read(output/'cad-slope/mesh.json')==read(cad/'mesh.json')
        write(output/'cad-slope/cad-input-link.json',dict(passed=True,mesh_contents_equal=True,
            cad_mesh_sha256=sha(cad/'mesh.json'),analysis_mesh_sha256=sha(output/'cad-slope/mesh.json'),
            cad_ground_sha256=sha(cad/'ground.json'),analysis_ground_sha256=sha(output/'cad-slope/ground.json'),
            scope='Exact captured CAD mesh contents. Analysis explicitly changes dry-reference unit weight to 16 and elastic-reference cohesion to 500; geometry/body IDs are unchanged. Flow and ONE_WAY settings are recorded separately.'))
        first=next(i for i,f in enumerate(run['frames']) if f['phase']=='transient');restarts+=restart(solver,output/'cad-slope',run,[1,first-1,first])
    report=dict(schema='soil-fem-one-way-cli-verification/1',passed=True,solver_sha256=sha(solver),restarts=restarts,cad_slope_exercised=bool(cad),profiles=list(profiles),genuine_mechanical_exhaustion_preserves_both_states=True,plastic_rejection_and_clean_accepted_times_exact=True,invalid_input_no_result=True,existing_output_preserved=True,
        maximum_independent_mechanical_residual=max(v['maximum_relative_mechanical_force_residual'] for v in checks.values()),elastic_comparisons=references,derived_point_omission_preserves_exact_primary_state=True,
        scope='Fixed active-domain frozen-porosity ONE_WAY MC/MCC. Initial equilibrium, separate fluid terms, physical-time rollback/restart, independent elastic mechanics given water heads. No excavation/BIRTH water transfer, wet JOINT or independent general 3D plastic solution is inferred.')
    write(output/'verification.json',report);print('ONE_WAY accepted:',len(restarts),'fresh-process restarts',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--cad',type=Path);a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve(),a.cad.resolve() if a.cad else None)
