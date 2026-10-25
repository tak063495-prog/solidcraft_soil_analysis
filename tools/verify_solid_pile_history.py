import hashlib
import json
from pathlib import Path
import subprocess
import sys

import argparse
parser=argparse.ArgumentParser();parser.add_argument('run',type=Path)
args=parser.parse_args()
folder=args.run.resolve()
root=Path(__file__).resolve().parents[1]
output=folder/'history';output.mkdir()
solver=folder/'solver/soil-fem-3d.exe'
mesh=folder/'coarse/mesh.json';ground=folder/'ground.json';case=folder/'stages.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
baseline=read(folder/'coarse/stage-results.json')
solver_digest=json.loads((folder/'solver-binding.json').read_text())['solver_sha256']
if hashlib.sha256(solver.read_bytes()).hexdigest()!=solver_digest or baseline['provenance']['solver_binary_sha256']!=solver_digest:
    raise ValueError('pile restart solver identity differs from frozen binding')
fingerprints={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (mesh,ground,case,solver,Path(__file__))}
def stable(x):
    if isinstance(x,dict):return {k:stable(v) for k,v in x.items() if k!='elapsed_seconds'}
    if isinstance(x,list):return [stable(v) for v in x]
    return x
geo=sum(a['accepted'] for a in baseline['geostatic']['attempts'])
loaded=sum(a['accepted'] and a['stage']=='load_pile' for a in baseline['attempts'])
reports=[]
for count in (1,geo+1,geo+loaded+1):
    checkpoint=output/f'checkpoint-{count}.json';result=output/f'resumed-{count}.json'
    for command,last,log in (
        ('checkpoint-stages',[str(count),str(checkpoint)],output/f'checkpoint-{count}.log'),
        ('resume-stages',[str(checkpoint),str(result)],output/f'resume-{count}.log'),
    ):
        subprocess.run([sys.executable,'-X','utf8',str(root/'tools/run_logged.py'),str(log),str(solver),command,
            str(mesh),str(ground),str(case),*last],cwd=root,check=True)
    saved=read(checkpoint);actual=read(result)
    assert saved['run']['failure_reason']=='checkpoint_requested'
    assert actual['complete'] and stable(actual)==stable(baseline)
    reports.append(dict(stop_after_accepted_increments=count,
        initial_progress=saved['run']['geostatic']['progress'],
        stage_index=saved['run']['last_committed']['stage_index'],stage_progress=saved['run']['last_committed']['stage_progress'],
        full_non_time_result_exact=True))
record=dict(passed=True,scope='finite solid single pile, gravity/load/unload; all non-time result fields and trial histories',
    checkpoints=reports,solver_sha256=solver_digest,input_and_source_sha256=fingerprints,
    spatial_convergence_accepted=False,goal_complete=False,p6_complete=False)
if fingerprints!={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in fingerprints}:
    raise ValueError('pile restart inputs or sources changed during verification')
with (output/'verification.json').open('x') as stream:json.dump(record,stream,indent=2)
print(json.dumps(record),flush=True)
