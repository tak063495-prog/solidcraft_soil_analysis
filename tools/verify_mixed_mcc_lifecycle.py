"""Separate-process checkpoints and converged-trial rejection for mixed MCC."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):
    with p.open('x',encoding='utf-8') as f:json.dump(d,f,indent=2)


def no_time(value):
    if isinstance(value,dict):return {k:no_time(v) for k,v in value.items() if k!='elapsed_seconds'}
    if isinstance(value,list):return [no_time(v) for v in value]
    return value


def exact(a,b):
    return json.dumps(no_time(a),sort_keys=True,allow_nan=False)==json.dumps(no_time(b),sort_keys=True,allow_nan=False)


def verify(base,out):
    root=Path(__file__).resolve().parents[1]
    solver=base/'solver/soil-fem-3d.exe';data=base/'n4'
    mesh=data/'mesh.json';ground=base/'ground.json';case=data/'stages-finest.json'
    native=read(data/'native-finest.json');assert native['complete']
    assert native['case']['schema']=='soil-fem-mixed-mcc-stages/1'
    assert sum(a['accepted'] for a in native['attempts'])==192
    inputs=[solver,mesh,ground,case,Path(__file__),root/'tools/run_logged.py']
    fingerprints={str(p):sha(p) for p in inputs}
    out.mkdir(exist_ok=False)
    def execute(tag,args,success=True):
        command=[sys.executable,'-X','utf8',str(root/'tools/run_logged.py'),str(out/(tag+'.log')),str(solver),*map(str,args)]
        process=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8')
        if (process.returncode==0)!=success:raise AssertionError(tag+': '+process.stdout[-2000:])
        print(tag+(' completed' if success else ' refused as expected'),flush=True)
    restarts=[]
    for count in (32,96,160):
        checkpoint=out/f'checkpoint-{count}.json';result=out/f'resumed-{count}.json'
        execute(f'checkpoint-{count}',['checkpoint-stages',mesh,ground,case,count,checkpoint])
        execute(f'resume-{count}',['resume-stages',mesh,ground,case,checkpoint,result])
        restored=read(result)
        assert restored['complete'] and exact(native,restored),count
        restarts.append(dict(accepted_increments=count,full_result_except_elapsed_seconds_exact=True))
    control=copy.deepcopy(read(case));control['settings']['initial_increment']=1/128
    injected=copy.deepcopy(read(case));injected['verification_reject_once_stages']=['load','unload','cut']
    write(out/'control-case.json',control);write(out/'injected-case.json',injected)
    execute('control',['solve-stages',mesh,ground,out/'control-case.json',out/'control.json'])
    execute('injected',['solve-stages',mesh,ground,out/'injected-case.json',out/'injected.json'])
    a,b=read(out/'control.json'),read(out/'injected.json')
    assert a['complete'] and b['complete']
    assert exact(a['frames'],b['frames']) and exact(a['last_committed'],b['last_committed'])
    rejected=[r for r in b['attempts'] if not r['accepted']]
    assert len(rejected)==3 and all(r['committed_state_unchanged_on_rejection'] for r in rejected)
    assert rejected[0]['rejected_trial_plastic_changed']
    assert rejected[-1]['rejected_trial_active_set_changed']
    changed=copy.deepcopy(read(case));changed['stages'][0]['target_pressure_kpa']+=1.
    write(out/'changed-load-case.json',changed)
    execute('wrong-load-refused',['resume-stages',mesh,ground,out/'changed-load-case.json',out/'checkpoint-32.json',out/'must-not-exist.json'],False)
    assert not (out/'must-not-exist.json').exists()
    assert fingerprints=={str(p):sha(p) for p in inputs}
    report=dict(schema='soil-fem-mixed-mcc-lifecycle-verification/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Dry two-layer mixed-MCC force cycle/full-layer removal; all primary pressure, cumulative volume and MCC histories',
        restarts=restarts,rejected_converged_trials=3,all_control_frames_and_states_exact=True,
        wrong_load_checkpoint_refused=True,input_and_source_sha256=fingerprints,
        evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()})
    write(out/'verification.json',report)
    print('Mixed-MCC lifecycle verified: 3 process restarts, 3 rejected converged trials, full histories exact.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('base',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();verify(a.base.resolve(),a.output.resolve())
