"""Reject converged load/reversal/cut trials and restart their committed histories.

The full and half-step paths have identical initial equilibrium. Only the trial
size and explicit rejection hooks differ; no production material update changes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def no_elapsed(value):
    if isinstance(value,dict):return {k:no_elapsed(v) for k,v in value.items() if k!='elapsed_seconds'}
    if isinstance(value,list):return [no_elapsed(v) for v in value]
    return value


def verify(mesh,ground,case_path,out):
    root=Path(__file__).resolve().parents[1];solver=root/'target/release/soil-fem-3d.exe'
    paths=[mesh,ground,case_path,solver,Path(__file__),root/'tools/run_logged.py']
    fingerprints={str(p):digest(p) for p in paths}
    original=read(case_path)
    if [s['kind'] for s in original['stages']]!=['boundary','boundary','excavate']:
        raise ValueError('this verification requires load, unload, excavation')
    if original['stages'][0]['target_factor']!=1 or original['stages'][1]['target_factor']!=0:
        raise ValueError('expected a complete boundary load/unload cycle')
    out.mkdir(exist_ok=False)
    def write(name,value):
        path=out/name
        with path.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2)
        return path
    def run(name,args):
        subprocess.run([sys.executable,str(root/'tools/run_logged.py'),str(out/(name+'.log')),
                        str(solver),*map(str,args)],cwd=root,check=True)
    variants={}
    for name,increment in (('rollback',1.),('control',.5)):
        case=json.loads(json.dumps(original))
        case['verification_reject_once_stages']=[s['name'] for s in case['stages']] if name=='rollback' else []
        for stage in case['stages']:stage['initial_increment']=increment
        variants[name]=write(name+'-case.json',case)
        run(name,['solve-stages',mesh,ground,variants[name],out/(name+'.json')])
    expected=read(out/'rollback.json');control=read(out/'control.json')
    assert expected['complete'] and control['complete']
    for key in ('frames','last_committed','reference_displacement_m'):
        assert expected[key]==control[key],key
    assert no_elapsed(expected['geostatic'])==no_elapsed(control['geostatic'])
    rejected=[a for a in expected['attempts'] if not a['accepted']]
    assert [a['stage'] for a in rejected]==[s['name'] for s in original['stages']]
    for attempt in rejected:
        assert attempt['reason']=='verification_injected_rejection_after_convergence'
        assert attempt['progress_from']==0 and attempt['progress_to']==1
        assert attempt['committed_state_unchanged_on_rejection']
        assert attempt['rejected_trial_material_changed'] and attempt['rejected_trial_load_changed']
    assert rejected[0]['rejected_trial_plastic_changed']
    assert rejected[-1]['rejected_trial_active_set_changed']
    accepted=[a for a in expected['attempts'] if a['accepted']]
    assert len(accepted)==6 and all(a['progress_to'] in (.5,1.) for a in accepted)
    offset=sum(a['accepted'] for a in expected['geostatic']['attempts'])
    reports=[]
    for index,attempt in enumerate(accepted[:-1],1):
        count=offset+index;name=f'restart-{count}'
        checkpoint=out/(name+'-checkpoint.json');result=out/(name+'-result.json')
        run(name+'-save',['checkpoint-stages',mesh,ground,variants['rollback'],count,checkpoint])
        saved=read(checkpoint);assert not saved['run']['complete']
        run(name+'-resume',['resume-stages',mesh,ground,variants['rollback'],checkpoint,result])
        restored=read(result);assert restored['complete']
        for key in ('frames','last_committed','reference_displacement_m'):
            assert restored[key]==expected[key],(count,key)
        for key in ('geostatic','attempts'):
            assert no_elapsed(restored[key])==no_elapsed(expected[key]),(count,key)
        reports.append(dict(accepted_steps=count,after_stage=attempt['stage'],progress_to=attempt['progress_to'],
                            all_history_frames_reference_and_nontiming_diagnostics_exact=True))
    assert fingerprints=={str(p):digest(p) for p in paths},'source/input changed during verification'
    files={p.name:digest(p) for p in out.iterdir() if p.is_file()}
    record=dict(passed=True,goal_complete=False,p6_complete=False,
                scope='Actual Tet10 slope: plastic boundary load, unload, local excavation; exact rollback and process restarts',
                all_frames_history_and_reference_exact=True,plastic_reload_rejection=True,
                converged_rejections=[{k:a[k] for k in ('stage','progress_from','progress_to','reason',
                    'committed_state_unchanged_on_rejection','rejected_trial_plastic_changed','rejected_trial_active_set_changed')} for a in rejected],
                separate_process_restarts=reports,excluded_from_diagnostic_equality=['elapsed_seconds'],
                input_and_source_sha256=fingerprints,evidence_sha256=files)
    write('verification.json',record)
    print(json.dumps({k:record[k] for k in ('passed','scope','plastic_reload_rejection','separate_process_restarts')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('ground',type=Path)
    p.add_argument('case',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    verify(*(v.resolve() for v in (a.mesh,a.ground,a.case,a.output)))
