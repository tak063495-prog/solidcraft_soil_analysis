"""Exact cutback replay of IDH shear/cap/void history on the accepted path."""
import argparse
import hashlib
import json
from pathlib import Path
from verify_checkpoint import invoke


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(solver,run,output):
    output.mkdir(parents=True,exist_ok=False);baseline=read(run/'stage-results.json')
    assert baseline['complete'] and baseline['schema']=='soil-fem-stage-result/10' and baseline.get('geostatic') is None
    assert baseline['provenance']['solver_binary_sha256']==sha(solver)
    rejected=[a for a in baseline['attempts'] if not a['accepted']]
    assert rejected and all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    assert any(a['stage']=='preload' and a['rejected_trial_plastic_changed'] for a in rejected)
    assert any(a['rejected_trial_active_set_changed'] and a['rejected_trial_load_changed'] for a in rejected)
    case=read(run/'stages.json');case['verification_reject_once_stages']=[]
    for stage in case['stages']:
        stage['initial_increment']=next(a['progress_to'] for a in baseline['attempts'] if a['accepted'] and a['stage']==stage['name'])
    cp=output/'clean-stages.json'
    with cp.open('x',encoding='utf-8') as f:json.dump(case,f,indent=2)
    rp=output/'clean-result.json'
    invoke([str(solver),'solve-stages',str(run/'ground/mesh.json'),str(run/'ground/ground.json'),str(cp),str(rp)],output/'clean.log')
    clean=read(rp);assert clean['complete']
    path=lambda d:[(a['stage'],a['progress_to']) for a in d['attempts'] if a['accepted']]
    assert path(clean)==path(baseline)
    for field in ('last_committed','frames','excavation_audits'):
        assert clean[field]==baseline[field],field
    assert clean.get('reference_displacement_m')==baseline.get('reference_displacement_m')
    report={'schema':'soil-fem-hardening-history-verification/1','passed':True,'solver_sha256':sha(solver),
        'exact_accepted_path_states_frames_and_audits':True,'rejected_trial_plastic_history_changed':True,
        'all_rejected_trials_preserved_committed_state':True,
        'excavation_rejected_trial_plastic_changed':any(a['rejected_trial_active_set_changed'] and a['rejected_trial_plastic_changed'] for a in rejected),
        'sha256':{p.name:sha(p) for p in output.iterdir() if p.is_file()}}
    with (output/'verification.json').open('x') as f:json.dump(report,f,indent=2)
    print('IDH rejected trial and clean accepted path match exactly:',run)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.run.resolve(),a.output_dir.resolve())
