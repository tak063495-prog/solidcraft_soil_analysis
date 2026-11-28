"""Small input-only crack-tip acceptance; preserves all raw runs and failures.

This is deliberately separate from historical JOINT acceptance tolerances.
Every new elastic comparison uses relative 1e-7 and force reassembly 1e-10.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import fem_reference as ref
from joint_reference import System
from joint_tip_fixture import fixture
from verify_joint_coupling import reference
from verify_high_order import prepare, internal, physical

ROOT = Path(__file__).resolve().parents[1]


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n',encoding='utf-8')


def invoke(command, log, success=True):
    with log.open('x',encoding='utf-8') as stream:
        stream.write(repr([str(p) for p in command])+'\n');stream.flush()
        result=subprocess.run([str(p) for p in command],stdout=stream,stderr=subprocess.STDOUT)
        stream.write(f'\nexit_code={result.returncode}\n')
    assert (result.returncode == 0) == success, (result.returncode, str(log))


def without_time(value):
    if isinstance(value,dict):return {k:without_time(v) for k,v in value.items() if k!='elapsed_seconds'}
    if isinstance(value,list):return [without_time(v) for v in value]
    return value


def solve(solver, folder, mesh, ground, case):
    folder.mkdir(parents=True,exist_ok=False)
    for name,data in [('mesh',mesh),('ground',ground),('stages',case)]:write(folder/(name+'.json'),data)
    inputs=[folder/'mesh.json',folder/'ground.json',folder/'stages.json']
    invoke([solver,'solve-stages',*inputs,folder/'result.json'],folder/'solve.log')
    run=json.loads((folder/'result.json').read_text())
    assert run['complete'],run.get('failure_reason')
    return run, inputs


def reassembly(mesh,ground,case,run):
    system=System(mesh,case['joint_system'],run['joint_integration_points'])
    q,offset=prepare(mesh,run['integration_points']);rows=[]
    linear=[r['true_relative_residual'] for a in run['geostatic']['attempts']+run['attempts'] for r in a.get('linear_solves',[])]
    assert linear and max(linear)<1e-10
    for frame in run['frames']:
        state=frame['state'];active=np.asarray(state['active_elements'])
        bulk,_=internal(mesh,ground,state['material_states'],active,q,offset,state.get('element_material_ids'),state.get('pore_pressure_kpa'))
        force=bulk+system.audit(state,frame['joint_results'])
        load,_=physical(mesh,ground,active,state['surface_pressure_kpa'],case['pressure_on_original_surface_only'],state.get('element_material_ids'),system)
        load+=state.get('boundary_factor',0)*ref.conditions_load(mesh,case['boundary'],active)
        applied=np.asarray(state['applied_nodal_load_kn']).reshape(-1,3)
        load_error=float(np.max(np.abs(applied-load)))
        assert load_error<1e-10,load_error
        basis,constraint_error=ref.staged_constraint_space(mesh,case,state)
        residual=np.concatenate([b.T@(force-applied)[n] for n,b in basis.items()])
        relative=float(np.linalg.norm(residual)/max(np.linalg.norm(applied),1.))
        assert relative<1e-10,(frame['name'],relative)
        assert constraint_error<1e-10 and frame['relative_force_residual']<1e-8
        support=system.rigid_support(case,state,frame['joint_support'])
        assert support['unrestrained_modes']==0
        rows.append(dict(frame=frame['name'],reassembled_relative_residual=relative,
                         load_max_abs_kn=load_error,constraint_error_m=float(constraint_error)))
    return rows


def compare_reference(mesh,ground,case,run):
    frames,residual=reference(mesh,ground,case,run,1)
    assert residual<1e-10
    materials={m['id']:m['mohr_coulomb'] for m in ground['materials']}
    body={r['id']:r['body'] for r in mesh['regions']}
    rows=[]
    for frame,expected in zip(run['frames'],frames):
        state=frame['state'];actual=np.asarray(state['displacement_m']).ravel()
        ue=float(np.linalg.norm(actual-expected['u'])/max(np.linalg.norm(expected['u']),1e-30))
        actual_q=np.asarray([s['traction_kpa'] for s in state['joints']['points']])
        qe=float(np.linalg.norm(actual_q-expected['traction'])/max(np.linalg.norm(expected['traction']),1.))
        history_errors={}
        for key,reference_key in [('plastic_slip_m','plastic'),('accumulated_plastic_slip_m','accumulated')]:
            actual_history=np.asarray([s[key] for s in state['joints']['points']]);reference_history=expected[reference_key]
            # An exactly elastic history is compared absolutely; nonzero slip
            # uses the same relative 1e-7 gate as the other independent fields.
            difference=float(np.linalg.norm(actual_history-reference_history));scale=float(np.linalg.norm(reference_history))
            error=difference/scale if scale else difference
            assert error<(1e-7 if scale else 1e-12),(frame['name'],key,error)
            history_errors[key]=error
        stresses=[]
        for e in range(len(mesh['tetrahedra'])):
            op=ref.volume(mesh,e);mat=materials[mesh['material_ids'][e]]
            seed=np.asarray(ground['initialization']['geostatic']['seed']['states'][body[mesh['region_ids'][e]]]['effective_stress_compression_positive_kpa'])
            ids=(3*op['nodes'][:,None]+np.arange(3)).ravel()
            stresses.extend(seed-np.einsum('ab,qbi,i->qa',ref.elastic_tensor(mat['young_kpa'],mat['poisson']),op['b'],expected['u'][ids]))
        actual_s=np.asarray([s['effective_stress_compression_positive_kpa'] for s in state['material_states']])
        se=float(np.linalg.norm(actual_s-stresses)/max(np.linalg.norm(stresses),1.))
        assert all(s['equivalent_plastic_strain']==0 for s in state['material_states'])
        assert max(ue,qe,se)<1e-7,(frame['name'],ue,qe,se)
        rows.append(dict(frame=frame['name'],all_node_displacement_relative_l2=ue,
                         all_volume_gp_stress_relative_l2=se,all_contact_gp_traction_relative_l2=qe,all_contact_gp_history_errors=history_errors))
    return dict(reference='input-only GL125 elastic 3D + independent Tri6/Coulomb Newton; one local substep in both implementations',
                maximum_independent_residual=residual,frames=rows,passed=True)


def restart_and_rejection(solver,folder,run,inputs):
    accepted=0;points=[];seen=set()
    # Count the initialization increments too: checkpoint CLI counts both.
    for a in run['geostatic']['attempts']:
        if a['accepted']:accepted+=1
    for a in run['attempts']:
        if not a['accepted']:continue
        accepted+=1
        if a['stage'] in seen or a['progress_to']>=1:continue
        seen.add(a['stage']);cp=folder/(a['stage']+'-checkpoint.json');rp=folder/(a['stage']+'-resumed.json')
        invoke([solver,'checkpoint-stages',*inputs,accepted,cp],folder/(a['stage']+'-checkpoint.log'))
        saved=json.loads(cp.read_text());assert saved['run']['last_committed']['stage_progress']<1
        invoke([solver,'resume-stages',*inputs,cp,rp],folder/(a['stage']+'-resume.log'))
        resumed=json.loads(rp.read_text());assert without_time(resumed)==without_time(run)
        points.append(dict(stage=a['stage'],accepted_increment=accepted,exact_except_elapsed_seconds=True))
    assert seen=={s['name'] for s in run['case']['stages']}
    rejected=[a for a in run['attempts'] if not a['accepted']]
    assert all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    injected=[a for a in rejected if a['reason']=='verification_injected_rejection_after_convergence']
    selected=run['case']['verification_reject_once_stages'];assert {a['stage'] for a in injected}==set(selected)
    case=json.loads(inputs[2].read_text());case['verification_reject_once_stages']=[]
    for stage in case['stages']:
        if stage['name'] in selected:stage['initial_increment']=stage.get('initial_increment',case['settings']['initial_increment'])/2
    write(folder/'clean-stages.json',case)
    invoke([solver,'solve-stages',*inputs[:2],folder/'clean-stages.json',folder/'clean.json'],folder/'clean.log')
    clean=json.loads((folder/'clean.json').read_text());assert clean['complete']
    assert clean['last_committed']==run['last_committed'] and clean['frames']==run['frames']
    return dict(restarts=points,intended_rejections=len(injected),clean_accepted_path_exact=True)


def main(solver,output):
    output.mkdir(parents=True,exist_ok=False)
    _,_,mesh,cfg=fixture()
    ground=json.loads((ROOT/'examples/ground-joint-crack-tip.json').read_text())
    template=json.loads((ROOT/'examples/stages-joint-crack-tip.json').read_text())
    cases={}
    for name,vector,factors in [('compression',[1e-5,0.,-1e-4],[1.]),('opening',[0.,0.,.02],[1.]),('reverse_shear',[.004,.001,0.],[1.,-1.])]:
        case=copy.deepcopy(template)
        for condition,v in zip(case['boundary']['constraints'][-3:],vector):condition['value_m']=v
        case['stages']=[dict(kind='boundary',name=name+str(i),target_factor=f,initial_increment=.5) for i,f in enumerate(factors)]
        case['verification_reject_once_stages']=[]
        folder=output/name;run,_=solve(solver,folder,mesh,ground,case)
        status=[s['status'] for s in run['frames'][-1]['state']['joints']['points']]
        expected={'compression':'stick','opening':'open','reverse_shear':'slip'}[name]
        assert status and all(s==expected for s in status)
        cases[name]=dict(assembly=reassembly(mesh,ground,run['case'],run),independent_3d=compare_reference(mesh,ground,run['case'],run),
                         final_contact_status=expected,contact_points=len(status),history=restart_and_rejection(solver,folder,run,[folder/'mesh.json',folder/'ground.json',folder/'stages.json']))
    case=copy.deepcopy(template);case['stages']=[dict(kind='excavate',name='remove_side',bodies=['UpperClay'],remaining_surface_pressure_kpa=0.,initial_increment=.5)]
    case['verification_reject_once_stages']=['remove_side']
    folder=output/'construction';run,inputs=solve(solver,folder,mesh,ground,case)
    construction=dict(assembly=reassembly(mesh,ground,run['case'],run),history=restart_and_rejection(solver,folder,run,inputs))
    assert all(not p['active'] for p in run['frames'][-1]['joint_results']['points'])
    assert run['frames'][0]['state']['joints']==run['frames'][-1]['state']['joints']
    assert not any(run['last_committed']['excavation_release_load_kn'])
    report=dict(schema='soil-fem-joint-crack-tip-acceptance/1',passed=True,dry_cases=cases,construction=construction,
                solver_sha256=hashlib.sha256(solver.read_bytes()).hexdigest())
    write(output/'measurements.json',report)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.solver.resolve(),a.output.resolve())
