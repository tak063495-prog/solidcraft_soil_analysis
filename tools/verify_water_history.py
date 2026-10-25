"""Exact replay after rejected hydrostatic trials, including old plasticity."""
import argparse
import json
from pathlib import Path
import hashlib
from verify_checkpoint import invoke


def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(solver,run,output,require_plastic=False):
    output.mkdir(parents=True,exist_ok=False);data=read(run/'stage-results.json');assert data['complete'] and data['schema'] in ('soil-fem-stage-result/11','soil-fem-stage-result/12')
    assert data['provenance']['solver_binary_sha256']==sha(solver)
    rejected=[a for a in data['attempts'] if not a['accepted']]
    assert rejected and all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    assert any(a.get('rejected_trial_water_changed') for a in rejected)
    if require_plastic:assert any(a['rejected_trial_plastic_changed'] for a in rejected)
    case=read(run/'stages.json');case['verification_reject_once_stages']=[]
    for stage in case['stages']:stage['initial_increment']=next(a['progress_to'] for a in data['attempts'] if a['accepted'] and a['stage']==stage['name'])
    cp=output/'clean-stages.json';rp=output/'clean-result.json'
    with cp.open('x') as f:json.dump(case,f,indent=2)
    invoke([str(solver),'solve-stages',str(run/'ground/mesh.json'),str(run/'ground/ground.json'),str(cp),str(rp)],output/'clean.log')
    clean=read(rp);assert clean['complete']
    path=lambda d:[(a['stage'],a['progress_to']) for a in d['attempts'] if a['accepted']]
    assert path(clean)==path(data)
    for field in ('frames','last_committed','excavation_audits'):assert clean[field]==data[field],field
    def stable_flows(rows):
        rows=json.loads(json.dumps(rows))
        for audit in rows:
            solution=audit['solution']
            for solve in ([solution['linear']] if solution['linear'] is not None else [])+solution.get('refinement_solves',[]):solve['elapsed_seconds']=0.
            for balance in solution['balances']:
                for solve in balance['linear_solves']:solve['elapsed_seconds']=0.
        return rows
    assert stable_flows(clean.get('flow_audits',[]))==stable_flows(data.get('flow_audits',[]))
    assert clean.get('reference_displacement_m')==data.get('reference_displacement_m')
    report={'schema':'soil-fem-water-history/1','passed':True,'solver_sha256':sha(solver),'exact_accepted_path_and_all_states':True,
        'rejected_water_trial_restored':True,'rejected_plastic_trial_exercised':any(a['rejected_trial_plastic_changed'] for a in rejected),
        'sha256':{p.name:sha(p) for p in output.iterdir() if p.is_file()}}
    with (output/'verification.json').open('x') as f:json.dump(report,f,indent=2)
    print('Staged water, plasticity and removed-element histories match the clean path exactly:',run)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--require-plastic',action='store_true')
    a=p.parse_args();verify(a.solver.resolve(),a.run.resolve(),a.output_dir.resolve(),a.require_plastic)
