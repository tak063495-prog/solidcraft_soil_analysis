"""Compare deliberate structural cutbacks with a clean accepted path."""
import argparse
import hashlib
import json
from pathlib import Path
from verify_checkpoint import invoke


def verify(solver,run,output):
    output.mkdir(parents=True,exist_ok=False)
    baseline=json.loads((run/'stage-results.json').read_text());assert baseline['complete']
    case=json.loads((run/'stages.json').read_text());selected=case['verification_reject_once_stages']
    rejects=[a for a in baseline['attempts'] if not a['accepted']]
    assert len(rejects)==len(selected)
    assert all(a['committed_state_unchanged_on_rejection'] and a['rejected_trial_structure_changed'] and a['rejected_trial_material_changed'] for a in rejects)
    mcc=any('mcc' in m for m in baseline['ground']['materials'])
    if mcc:assert any(a['rejected_trial_plastic_changed'] for a in rejects)
    for stage in case['stages']:
        if stage['name'] in selected:
            stage['initial_increment']=(stage.get('initial_increment') or case['settings']['initial_increment'])/2
    case['verification_reject_once_stages']=[]
    cp=output/'clean-stages.json';cp.write_text(json.dumps(case,indent=2));result=output/'clean-result.json'
    invoke([str(solver),'solve-stages',str(run/'ground/mesh.json'),str(run/'ground/ground.json'),str(cp),str(result)],output/'clean.log')
    clean=json.loads(result.read_text());assert clean['complete']
    assert clean['last_committed']==baseline['last_committed']
    assert clean['frames']==baseline['frames']
    report={'schema':'soil-fem-beam-history-verification/1','solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),
            'exact_final_and_all_frames':True,'rejected_structure_and_material_changed':True,
            'rejected_plastic_changed':mcc,'rejected_stages':[a['stage'] for a in rejects],
            'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (run/'stage-results.json',cp,result,Path(__file__))}}
    with (output/'verification.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.run.resolve(),a.output_dir.resolve())
