"""Actual CLI numerical failure and policy rejection must not become Fs bounds."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
from verify_checkpoint import invoke
from export_replays import capture


def verify(solver,source,output):
    output.mkdir(parents=True,exist_ok=False);mesh=source/'ground/mesh.json'
    ground=json.loads((source/'ground/ground.json').read_text());case=json.loads((source/'stages.json').read_text())
    records=[]
    for kind in ('numerical','policy'):
        folder=output/kind;folder.mkdir();g=copy.deepcopy(ground);c=copy.deepcopy(case)
        spec={'schema':'soil-fem-replay/1','material_ids':[m['id'] for m in g['materials']],'factors':[1.],
            'psi_policy':'keep','tension_policy':'keep'}
        if kind=='numerical':
            c['settings'].update(maximum_newton_iterations=1,initial_increment=1.,minimum_increment=1.)
            for m in g['materials']:m['mohr_coulomb'].update(cohesion_kpa=1.,friction_deg=10.,dilation_deg=5.)
        else:spec['factors']=[1.,100.]
        for name,value in (('ground.json',g),('stages.json',c),('specification.json',spec)):
            with (folder/name).open('x') as f:json.dump(value,f,indent=2)
        result=folder/'replay-results.json'
        invoke([str(solver),'replay-stages',str(mesh),str(folder/'ground.json'),str(folder/'stages.json'),str(folder/'specification.json'),str(result)],folder/'replay.log',1)
        native=json.loads(result.read_text());assert native['stability_classification']=='not_evaluated'
        if kind=='numerical':
            trial=native['trials'][0];assert trial['outcome']=='incomplete';run=trial['run'];assert not run['complete'] and not run['frames']
            assert run['geostatic']['seed']==run['last_committed']
            assert any('Newton iteration limit' in a['reason'] for a in run['geostatic']['attempts'])
            assert all(a['committed_state_unchanged_on_rejection'] for a in run['geostatic']['attempts'])
        else:
            assert native['trials'][0]['outcome']=='completed' and native['trials'][1]['outcome']=='rejected_input'
            assert native['trials'][1]['adjusted_ground'] is None and native['trials'][1]['run'] is None
        capture(mesh,result,folder/'trials',solver)
        package=json.loads((folder/'trials/verification.json').read_text())
        for trial,row in zip(native['trials'],package['trials']):
            assert row['accepted_export']==(trial['outcome']=='completed')
            if not row['accepted_export']:assert not (folder/'trials'/row['folder']/'results').exists()
        records.append({'case':kind,'outcomes':[t['outcome'] for t in native['trials']],
            'noncompleted_accepted_exports':False,'stability_or_fs_inferred':False})
    report={'schema':'soil-fem-replay-failure-verification/1','solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),'cases':records,
        'sha256':{str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in output.rglob('*.json')}}
    with (output/'verification.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.run.resolve(),a.output_dir.resolve())
