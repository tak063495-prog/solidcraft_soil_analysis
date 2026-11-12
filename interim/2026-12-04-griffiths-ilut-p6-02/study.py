"""Bounded fresh production solves, measured comparisons; no tolerance fitting."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import time
import numpy as np
from prepare import ROOT, prepare, read, write


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def bits(v):
    if isinstance(v, dict): return {k:bits(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v, list): return [bits(x) for x in v]
    if isinstance(v, float): return ('f64',struct.pack('>d',v).hex())
    return v


def command(args, log, budget=180):
    start=time.monotonic()
    with log.open('x',encoding='utf-8') as f:
        p=subprocess.Popen([str(x) for x in args],stdout=f,stderr=subprocess.STDOUT,cwd=ROOT)
        timeout=False
        try: code=p.wait(timeout=budget)
        except subprocess.TimeoutExpired:
            timeout=True; p.terminate()
            try: code=p.wait(timeout=10)
            except subprocess.TimeoutExpired: p.kill();code=p.wait()
    return dict(argv=[str(x) for x in args],pid=p.pid,exit_code=code,timed_out=timeout,
                elapsed_seconds=time.monotonic()-start,log_sha256=sha(log))


def metrics(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('full field shape/finite mismatch')
    d=a-b; norm=float(np.linalg.norm(a)); absolute=float(np.linalg.norm(d))
    return dict(shape=list(a.shape),components=a.size,maximum_absolute=float(abs(d).max(initial=0)),
                absolute_l2=absolute,relative_l2=absolute/norm if norm else (0. if not absolute else None),
                comparison_criterion=None,verdict='判定基準なし')


def fields(a,b):
    return dict(displacement_m=metrics(a['displacement_m'],b['displacement_m']),
        stress_kpa=metrics([s['effective_stress_compression_positive_kpa'] for s in a['material_states']],
                          [s['effective_stress_compression_positive_kpa'] for s in b['material_states']]))


def summary(run):
    rows=[]
    for stage,attempts in [('geostatic',run.get('geostatic',{}).get('attempts',[])),
                           *[(s['name'],[a for a in run['attempts'] if a['stage']==s['name']]) for s in run['case']['stages']]]:
        solves=[s for a in attempts for s in a['linear_solves']]
        accepted=[a for a in attempts if a['accepted']]
        rows.append(dict(stage=stage,attempts=len(attempts),accepted=len(accepted),
            newton_iterations=sum(a['newton_iterations'] or 0 for a in attempts),
            unreported_newton_attempts=sum(a['newton_iterations'] is None for a in attempts),
            maximum_newton_iterations=max((a['newton_iterations'] or 0 for a in attempts),default=0),
            gmres_iterations=sum(s['iterations'] for s in solves),
            maximum_gmres_iterations=max((s['iterations'] for s in solves),default=0),
            reported_linear_solves=len(solves),
            maximum_true_linear_residual=max((s['true_relative_residual'] for s in solves),default=0.),
            final_relative_force_residual=accepted[-1]['relative_force_residual'] if accepted else None,
            rejected_reasons=[a['reason'] for a in attempts if not a['accepted']]))
    return dict(complete=run['complete'],failure_reason=run['failure_reason'],stages=rows,
        force_gate=run['case']['settings']['relative_force_tolerance'],
        true_linear_gate=run['case']['settings']['relative_linear_tolerance'])


def compare(a,b):
    if len(a['frames'])!=len(b['frames']): raise ValueError('frame inventory differs')
    rows=[]
    for x,y in zip(a['frames'],b['frames'],strict=True):
        if x['name']!=y['name']: raise ValueError('frame identity differs')
        rows.append(dict(name=x['name'],**fields(x['state'],y['state'])))
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('output',type=Path)
    p.add_argument('--cli',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--strides',type=int,nargs='+',default=[5,2]);p.add_argument('--budget',type=int,default=180)
    a=p.parse_args(); out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    cli=a.cli.resolve();base=a.baseline.resolve();rows=[]
    for stride in a.strides:
        folder=out/f'stride{stride}';folder.mkdir()
        mesh,ground,case=prepare(folder/'input',stride)
        ilut=copy.deepcopy(case);ilut['settings']['linear_solver']['options']['ilut']=dict(drop_tolerance=1e-4,fill_limit=200)
        write(folder/'input/case-ilut.json',ilut)
        processes={};runs={}
        for label,binary,cp in [('base',base,'case.json'),('default',cli,'case.json'),('ilut',cli,'case-ilut.json')]:
            result=folder/(label+'.json')
            processes[label]=command([binary,'solve-stages',folder/'input/mesh.json',folder/'input/ground.json',
                folder/'input'/cp,result],folder/(label+'.log'),a.budget)
            if result.exists(): runs[label]=read(result)
        row=dict(stride=stride,elements=len(mesh['tetrahedra']),nodes=len(mesh['nodes']),processes=processes,
                 summaries={k:summary(r) for k,r in runs.items()},default_unchanged=None,field_differences=None)
        if 'base' in runs and 'default' in runs:
            before=copy.deepcopy(runs['base']);after=copy.deepcopy(runs['default'])
            assert before['provenance']['solver_binary_sha256']==sha(base)
            assert after['provenance']['solver_binary_sha256']==sha(cli)
            before['provenance']['solver_binary_sha256']=after['provenance']['solver_binary_sha256']
            row['default_unchanged']=bits(before)==bits(after)
        if 'default' in runs and 'ilut' in runs:
            row['field_differences']=compare(runs['default'],runs['ilut'])
        rows.append(row)
        write(folder/'comparison.json',row)
    write(out/'measurements.json',dict(passed=False,p6_complete=False,p602_complete=False,cases=rows,
        binaries_sha256=dict(base=sha(base),current=sha(cli)),full_field_criterion='判定基準なし',
        unmet=['元5%基準','正確な斜面Fs','実断面SRM','xy2']))


if __name__=='__main__': main()
