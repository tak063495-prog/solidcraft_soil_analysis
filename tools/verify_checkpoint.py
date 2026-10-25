"""Exercise fresh CLI processes, exact history recovery and rejected inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def invoke(command,log,expected=0):
    with log.open('x',encoding='utf-8') as stream:
        stream.write('Command: '+repr(command)+'\n');stream.flush()
        result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
        stream.write(f'\nExit code: {result.returncode}\n')
    assert result.returncode==expected,(result.returncode,log.read_text(encoding='utf-8'))


def verify(solver,ground_run,output):
    output.mkdir(parents=True,exist_ok=False)
    for name in ('mesh.json','ground.json'):shutil.copyfile(ground_run/name,output/name)
    root=Path(__file__).resolve().parents[1];case=json.loads((root/'examples/stages-local-sparse.json').read_text())
    case['stages']=[
        {'kind':'pressure','name':'preload','target_pressure_kpa':90.,'initial_increment':1.},
        {'kind':'pressure','name':'unload','target_pressure_kpa':30.},
        case['stages'][0]]
    case['verification_reject_once_stages']=['preload','excavate_pit']
    case_path=output/'stages.json';case_path.write_text(json.dumps(case,indent=2))
    inputs=[str(output/'mesh.json'),str(output/'ground.json'),str(case_path)]
    baseline_path=output/'continuous.json';checkpoint_path=output/'checkpoint.json';resumed_path=output/'resumed.json'
    invoke([str(solver),'solve-stages',*inputs,str(baseline_path)],output/'continuous.log')
    baseline=json.loads(baseline_path.read_text());assert baseline['complete']
    # A converged trial is deliberately rejected *after* changing its material
    # history. Compare with a fresh process following the accepted step sizes;
    # this catches a trial history accidentally reused after cutback.
    rejected_attempts=[a for a in baseline['attempts'] if not a['accepted']]
    assert len(rejected_attempts)==2
    assert all(a['committed_state_unchanged_on_rejection'] and a['rejected_trial_material_changed']
               and a['rejected_trial_load_changed'] for a in rejected_attempts)
    assert any(a['rejected_trial_active_set_changed'] for a in rejected_attempts)
    assert any(a['stage']=='preload' and a['rejected_trial_plastic_changed']
               for a in rejected_attempts),'preload rejection must exercise changed plastic history'
    clean_case=json.loads(json.dumps(case));clean_case['verification_reject_once_stages']=[]
    for stage in clean_case['stages']:
        if stage['name'] in case['verification_reject_once_stages']:
            stage['initial_increment']=stage.get('initial_increment',case['settings']['initial_increment'])/2
    clean_case_path=output/'clean-stages.json';clean_case_path.write_text(json.dumps(clean_case,indent=2))
    clean_path=output/'clean.json'
    invoke([str(solver),'solve-stages',*inputs[:2],str(clean_case_path),str(clean_path)],output/'clean.log')
    clean=json.loads(clean_path.read_text());assert clean['complete']
    assert baseline['last_committed']==clean['last_committed'],'trial history leaked through cutback'
    assert [f['state'] for f in baseline['frames']]==[f['state'] for f in clean['frames']]
    count=0
    for attempt in baseline['attempts']:
        if attempt['accepted']:
            count+=1
            if attempt['stage']=='excavate_pit':break
    invoke([str(solver),'checkpoint-stages',*inputs,str(count),str(checkpoint_path)],output/'checkpoint.log')
    saved=json.loads(checkpoint_path.read_text());state=saved['run']['last_committed']
    assert saved['schema']=='soil-fem-checkpoint/1'
    assert not saved['run']['complete'] and state['stage_index']==3 and 0<state['stage_progress']<1
    assert any(abs(x)>0 for x in state['excavation_release_load_kn'])
    assert any(any(abs(x)>1e-8 for x in p['plastic_strain_compression_positive']) for p in state['material_states'])
    assert saved['identity']['solver_binary_sha256']==hashlib.sha256(solver.read_bytes()).hexdigest()
    invoke([str(solver),'resume-stages',*inputs,str(checkpoint_path),str(resumed_path)],output/'resume.log')
    baseline=json.loads(baseline_path.read_text());resumed=json.loads(resumed_path.read_text())
    assert baseline['complete'] and resumed['complete']
    assert baseline['last_committed']==resumed['last_committed']
    assert baseline['frames']==resumed['frames']
    assert [(a['accepted'],a['progress_from'],a['progress_to']) for a in baseline['attempts']]==[
        (a['accepted'],a['progress_from'],a['progress_to']) for a in resumed['attempts']]
    # Check the other load-path phases in fresh processes as well. The existing
    # checkpoint above specifically exercises a nonzero excavation holding load.
    restarts=[{'stage':'excavate_pit','accepted_increment':count,
               'saved_progress':state['stage_progress'],'exact':True}]
    accepted=0;seen={'excavate_pit'}
    for attempt in baseline['attempts']:
        if not attempt['accepted']:continue
        accepted+=1
        name=attempt['stage']
        if name in seen:continue
        seen.add(name)
        cp=output/f'{name}-checkpoint.json';rp=output/f'{name}-resumed.json'
        invoke([str(solver),'checkpoint-stages',*inputs,str(accepted),str(cp)],output/f'{name}-checkpoint.log')
        snapshot=json.loads(cp.read_text());progress=snapshot['run']['last_committed']['stage_progress']
        assert not snapshot['run']['complete'] and 0<progress<1
        invoke([str(solver),'resume-stages',*inputs,str(cp),str(rp)],output/f'{name}-resume.log')
        restored=json.loads(rp.read_text());assert restored['complete']
        for key in ('last_committed','frames','excavation_audits'):
            assert restored[key]==baseline[key],(name,key)
        assert [(a['accepted'],a['progress_from'],a['progress_to']) for a in restored['attempts']]==[
            (a['accepted'],a['progress_from'],a['progress_to']) for a in baseline['attempts']]
        restarts.append({'stage':name,'accepted_increment':accepted,'saved_progress':progress,'exact':True})
    assert seen=={'preload','unload','excavate_pit'}
    # Element death archives every integration-point history, including pc/e
    # and plastic strain, rather than continuing to integrate removed soil.
    before_cut=next(f['state'] for f in baseline['frames'] if f['name']=='unload')
    removed_points=[i for i,p in enumerate(baseline['integration_points'])
                    if not baseline['last_committed']['active_elements'][p['element_index']]]
    assert removed_points
    assert all(before_cut['material_states'][i]==baseline['last_committed']['material_states'][i]
               for i in removed_points),'removed soil histories changed after excavation'
    # A modified mesh identity, even at the same sizes, must be refused.
    mesh=json.loads((output/'mesh.json').read_text());mesh['crs']='CHANGED'
    wrong=output/'changed-mesh.json';wrong.write_text(json.dumps(mesh))
    rejected=output/'rejected-output.json'
    invoke([str(solver),'resume-stages',str(wrong),*inputs[1:],str(checkpoint_path),str(rejected)],output/'changed-input.log',1)
    assert not rejected.exists()
    tampered=json.loads(checkpoint_path.read_text());tampered['run']['last_committed']['excavation_release_load_kn'][0]+=1.
    bad=output/'tampered-checkpoint.json';bad.write_text(json.dumps(tampered))
    invoke([str(solver),'resume-stages',*inputs,str(bad),str(rejected)],output/'tampered.log',1)
    assert not rejected.exists()
    # The original output must remain byte-identical on a repeated invocation.
    digest=hashlib.sha256(resumed_path.read_bytes()).hexdigest()
    invoke([str(solver),'resume-stages',*inputs,str(checkpoint_path),str(resumed_path)],output/'no-overwrite.log',1)
    assert hashlib.sha256(resumed_path.read_bytes()).hexdigest()==digest
    assert not list(output.glob('.*.tmp'))
    report={'schema':'soil-fem-checkpoint-verification/1','separate_process_resume':True,'exact_final_committed_state':True,
            'cutback_exact_vs_clean_accepted_path':True,'cutback_rejected_changed_material_load_and_active_set':True,
            'exact_frames':True,'exact_increment_path':True,'partial_excavation_progress':state['stage_progress'],
            'plastic_history_preserved':True,'changed_inputs_rejected':True,'corrupted_history_rejected':True,
            'rejected_trial_plastic_history_changed':True,'restarts':restarts,
            'removed_integration_point_histories_frozen':len(removed_points),
            'existing_output_unchanged':True,'temporary_files_remaining':False,
            'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),
            'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob('*.json')}}
    (output/'checkpoint-verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,default=Path('target/release/soil-fem-3d.exe'))
    p.add_argument('--ground-run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.ground_run.resolve(),a.output_dir.resolve())
