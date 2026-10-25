"""Separate-process restarts during equilibration, at its end and construction."""
import argparse
import hashlib
import json
from pathlib import Path
from verify_checkpoint import invoke
from export_stages import export


def stable_report(report):
    report=json.loads(json.dumps(report))
    for attempt in report['attempts']:
        for solve in attempt['linear_solves']:solve['elapsed_seconds']=0.
        for predictor in attempt.get('contact_predictors',[]):predictor['temporary_linear_solve']['elapsed_seconds']=0.
    return report


def verify(solver,run,output):
    output.mkdir(parents=True,exist_ok=False);rp=run/'stage-results.json';baseline=json.loads(rp.read_text())
    assert baseline['complete'] and baseline['schema']=='soil-fem-stage-result/4'
    assert baseline['provenance']['solver_binary_sha256']==hashlib.sha256(solver.read_bytes()).hexdigest()
    rejected=[a for a in baseline['attempts'] if not a['accepted']]
    assert all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    paths=[run/'ground/mesh.json',run/'ground/ground.json',run/'stages.json'];inputs=list(map(str,paths))
    count=sum(a['accepted'] for a in baseline['geostatic']['attempts']);reports=[]
    for stop in (1,count,count+1):
        cp=output/f'checkpoint-{stop}.json';resumed=output/f'resumed-{stop}.json'
        invoke([str(solver),'checkpoint-stages',*inputs,str(stop),str(cp)],output/f'checkpoint-{stop}.log')
        saved=json.loads(cp.read_text());assert saved['run']['failure_reason']=='checkpoint_requested'
        invoke([str(solver),'resume-stages',*inputs,str(cp),str(resumed)],output/f'resume-{stop}.log')
        actual=json.loads(resumed.read_text());assert actual['complete']
        assert actual['last_committed']==baseline['last_committed'] and actual['frames']==baseline['frames']
        assert actual['reference_displacement_m']==baseline['reference_displacement_m']
        assert stable_report(actual['geostatic'])==stable_report(baseline['geostatic'])
        reports.append({'stop_after_accepted_increments':stop,'geostatic_progress':saved['run']['geostatic']['progress'],
            'geostatic_complete':saved['run']['geostatic']['complete'],'exact_histories_frames_reference_and_increment_path':True})
    tampered=json.loads((output/'checkpoint-1.json').read_text());tampered['run']['geostatic']['initial_holding_load_kn'][0]+=1.
    cp=output/'tampered.json';cp.write_text(json.dumps(tampered));bad=output/'bad-result.json'
    invoke([str(solver),'resume-stages',*inputs,str(cp),str(bad)],output/'tampered.log',1);assert not bad.exists()
    ground=json.loads(paths[1].read_text());ground['initialization']['geostatic']['maximum_initial_engineering_strain']=1e-8
    gp=output/'too-large-initial-strain.json';gp.write_text(json.dumps(ground));failed=output/'failed-result.json'
    invoke([str(solver),'solve-stages',str(paths[0]),str(gp),str(paths[2]),str(failed)],output/'strain-limit.log',1)
    failure=json.loads(failed.read_text());assert not failure['complete'] and not failure['frames']
    assert failure['geostatic']['seed']==failure['geostatic']['last_committed']
    try:
        export(paths[0],failed,output/'failed-export');raise AssertionError('failed initialization exported')
    except ValueError as e:assert 'completed' in str(e)
    assert not (output/'failed-export').exists()
    report={'schema':'soil-fem-geostatic-checkpoint-verification/1','restarts':reports,
        'construction_rejected_trials_restored':len(rejected),
        'modified_holding_load_rejected':True,'strain_limit_preserves_seed_and_rejects_export':True,
        'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),
        'artifacts_sha256':{str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob('*.json')}}
    (output/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='artifacts_sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.run.resolve(),a.output_dir.resolve())
