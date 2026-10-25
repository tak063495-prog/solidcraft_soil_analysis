"""Rejected sliding histories must equal a fresh clean accepted load path."""
import argparse
import hashlib
import json
from pathlib import Path
from verify_checkpoint import invoke


def verify(solver,run,output):
    output.mkdir(parents=True,exist_ok=False)
    native=json.loads((run/'stage-results.json').read_text());assert native['complete']
    case=json.loads((run/'stages.json').read_text());selected=case['verification_reject_once_stages'];assert selected
    rejected=[a for a in native['attempts'] if not a['accepted']]
    injected=[a for a in rejected if a['reason']=='verification_injected_rejection_after_convergence']
    assert {a['stage'] for a in injected}==set(selected)
    assert all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    assert all(a.get('rejected_trial_structure_changed') and a['rejected_trial_material_changed'] for a in injected)
    assert any(a.get('rejected_trial_joint_slip_changed') for a in injected)
    case['verification_reject_once_stages']=[]
    for stage in case['stages']:
        if stage['name'] in selected:stage['initial_increment']=(stage.get('initial_increment') or case['settings']['initial_increment'])/2
    cp=output/'clean-stages.json';cp.write_text(json.dumps(case,indent=2),encoding='utf-8');rp=output/'clean-result.json'
    invoke([str(solver),'solve-stages',str(run/'ground/mesh.json'),str(run/'ground/ground.json'),str(cp),str(rp)],output/'clean.log')
    clean=json.loads(rp.read_text());assert clean['complete']
    assert native['last_committed']==clean['last_committed'] and native['frames']==clean['frames']
    report={'schema':'soil-fem-joint-history-verification/1','exact_final_and_all_frames':True,
        'rejected_trial_joint_slip_changed':True,'natural_rejections':len(rejected)-len(injected),
        'rejected_bulk_plastic_changed':any(a['rejected_trial_plastic_changed'] for a in injected),
        'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),
        'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (run/'stage-results.json',cp,rp,Path(__file__))}}
    with (output/'verification.json').open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.run.resolve(),a.output_dir.resolve())
