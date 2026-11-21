"""Tighter input-only stopping experiments, never new acceptance conditions.

Original nine-case inputs/defaults/gates are untouched. A resource timeout is
retained as incomplete and must not be classified as numerical convergence.
"""
import argparse
import copy
import json
import shutil
import subprocess
import time
from pathlib import Path
from common import ROOT,read,write,sha

def run(source,output,cli,budget):
    output.mkdir(parents=True,exist_ok=False)
    original=read(source/'stages.json');records=[]
    for name,changes in [('linear-only',{'relative_linear_tolerance':1e-12}),('force-and-linear',{'relative_linear_tolerance':1e-12,'relative_force_tolerance':1e-10})]:
        folder=output/name;folder.mkdir();case=copy.deepcopy(original)
        for key,value in changes.items():
            assert value<case['settings'][key];case['settings'][key]=value
        actual={k:(original['settings'][k],case['settings'][k]) for k in original['settings'] if original['settings'][k]!=case['settings'][k]}
        assert set(actual)==set(changes)
        for key in ('mesh','ground'):shutil.copy2(source/(key+'.json'),folder/(key+'.json'))
        write(folder/'stages.json',case)
        args=[str(cli),'solve-stages',*[str(folder/(s+'.json')) for s in ('mesh','ground','stages')],str(folder/'run.json')]
        started=time.monotonic();timeout=False
        with (folder/'process.log').open('x') as log:
            p=subprocess.Popen(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            try:rc=p.wait(timeout=budget)
            except subprocess.TimeoutExpired:
                timeout=True;p.terminate()
                try:rc=p.wait(timeout=20)
                except subprocess.TimeoutExpired:p.kill();rc=p.wait()
        result=dict(name=name,command=args,pid=p.pid,exit_code=rc,elapsed_seconds=time.monotonic()-started,
                    resource_budget_seconds=budget,resource_timeout=timeout,only_tighter_input_solver_settings=actual,production_defaults_changed=False,existing_acceptance_gates_changed=False,accepted_as_B02=False,
                    input_sha256={s:sha(folder/(s+'.json')) for s in ('mesh','ground','stages')},solver_sha256=sha(cli))
        if (folder/'run.json').exists():
            r=read(folder/'run.json');attempts=r['attempts']+r.get('geostatic',{}).get('attempts',[])
            result.update(complete=r['complete'],failure_reason=r['failure_reason'],frame_relative_force_residuals={f['name']:f['relative_force_residual'] for f in r['frames']},
                          maximum_accepted_force_residual=max((a['relative_force_residual'] or 0. for a in attempts if a['accepted']),default=None),
                          maximum_true_linear_residual=max((v['true_relative_residual'] for a in attempts for v in a.get('linear_solves',[])),default=None),
                          accepted_attempts=sum(a['accepted'] for a in attempts),rejected_attempts=sum(not a['accepted'] for a in attempts))
        else:result.update(complete=False,failure_reason='resource_timeout' if timeout else 'no_run')
        records.append(result);write(output/'stopping-controls.json',dict(passed=False,b02_complete=False,controls=records))
        print(name,rc,result['complete'],result.get('maximum_accepted_force_residual'),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('output',type=Path);p.add_argument('--cli',type=Path,required=True);p.add_argument('--budget',type=float,default=180.);a=p.parse_args();run(a.source.resolve(),a.output.resolve(),a.cli.resolve(),a.budget)
