"""C-01 global contact-switch cutback; inputs only, finite-step evidence, no extrapolation."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OLD=ROOT/'interim/2026-11-15-wet-joint-convergence'
sys.path.insert(0,str(OLD))
import study as prior
from study import read,write,sha,Writer,states,STAGES,FIELDS,fem
spec=importlib.util.spec_from_file_location('c01_global_reference',HERE/'reference.py')
reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
STEPS=(16,32,64,128,256)
BISECTIONS=4  # fixed before measuring; resolution only, never an acceptance tolerance
BASE='9338204990ebb529b811918dc173af7aa9572903'


def without_clock(v):
    if isinstance(v,dict): return {k:without_clock(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list): return [without_clock(x) for x in v]
    return v


def initialize(out,solver):
    out.mkdir(parents=True,exist_ok=False)
    for name,source in [('mesh-h0.json','joint-water-column-tet10-mesh.json'),('ground.json','ground-joint-water-mc.json'),('stages.json','stages-joint-water-cycle.json')]:
        value=read(ROOT/'examples'/source)
        if name=='ground.json': value=fem.normalize_ground(value)
        if name=='stages.json': value=fem.normalize_stage_case(value)
        write(out/name,value)
    binary=out/('probe.exe' if solver.suffix=='.exe' else 'probe');shutil.copy2(solver,binary)
    write(out/'identity.json',dict(starting_commit=BASE,solver_sha256=sha(binary),inputs={n:sha(out/n) for n in ('mesh-h0.json','ground.json','stages.json')},switch_bisections=BISECTIONS))


def audit(folder,mesh,ground,case,steps):
    spec=importlib.util.spec_from_file_location('c01_global_field_audit',HERE/'audit.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    result=module.audit(folder,mesh,ground,case,'cut',steps)
    coords=[(row['stage'],row['progress']) for row,_ in states(folder)]
    assert len(coords)==len(set(coords)) and coords==sorted(coords,key=lambda x:(STAGES.index(x[0]),x[1]))
    for stage in STAGES:
        assert {(stage,i/steps) for i in range(1,steps+1)} <= set(coords), 'prescribed base boundary missing'
    return result


def native(out,g,timeout,enabled=True):
    folder=out/f'g{g}';folder.mkdir(exist_ok=True)
    mesh=read(out/'mesh-h0.json');ground=read(out/'ground.json');case=read(out/'stages.json')
    cfg=prior.configuration(case,g,5e-5)
    if enabled:
        cfg['joint_system']['global_step_control']={'switch_bisections':BISECTIONS}
        cfg['settings']['minimum_increment']=min(cfg['settings']['minimum_increment'],1/(g*2**(BISECTIONS+1)))
    write(folder/'native-stages.json',cfg)
    binary=out/('probe.exe' if (out/'probe.exe').exists() else 'probe')
    cmd=[str(binary),str(out/'mesh-h0.json'),str(out/'ground.json'),str(folder/'native-stages.json'),str(folder/'native-raw')]
    started=time.monotonic()
    with (folder/'native.log').open('x') as log:
        r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
    assert r.returncode==0,('native failed',g,r.returncode)
    raw=folder/'native-raw';run=read(raw/'run.json');assert run['complete'] and len(run['frames'])==7
    phases=[('geostatic',run['geostatic']['attempts'])]+[(s['name'],[a for a in run['attempts'] if a['stage']==s['name']]) for s in cfg['stages']]
    summary=[]
    for name,attempts in phases:
        accepted=[a for a in attempts if a['accepted']]
        assert {i/g for i in range(1,g+1)} <= {a['progress_to'] for a in accepted}
        assert all(a['relative_force_residual'] <= cfg['settings']['relative_force_tolerance'] for a in accepted)
        assert all(a['committed_state_unchanged_on_rejection'] for a in attempts if not a['accepted'])
        assert all(x['true_relative_residual'] <= cfg['settings']['relative_linear_tolerance'] for a in attempts for x in a['linear_solves'])
        summary.append(dict(stage=name,accepted=len(accepted),rejected=len(attempts)-len(accepted),switch_cutbacks=sum(a['reason']=='joint_global_switch_cutback' for a in attempts),maximum_force_residual=max(a['relative_force_residual'] for a in accepted)))
    op=prior.operators(mesh,ground,case,'cut');writer=Writer(folder/'native')
    accepted_coordinates=[(name,a['progress_to']) for name,attempts in phases for a in attempts if a['accepted']]
    with (raw/'accepted.jsonl').open() as stream:
        for line,(name,progress) in itertools.zip_longest(stream,accepted_coordinates):
            assert line is not None and name is not None
            state=json.loads(line)
            writer(dict(stage=name,progress=progress,factor=state['boundary_factor'],status=[{'open':0,'stick':1,'slip':2}[p['status']] for p in state['joints']['points']]),prior.native_fields(state,op))
    writer.close()
    write(folder/'native-audit.json',dict(physical=audit(folder/'native',mesh,ground,case,g),phases=summary,elapsed_seconds=time.monotonic()-started,command=cmd))
    print('native complete',g,flush=True)


def independent(out,g,enabled=True):
    folder=out/f'g{g}';folder.mkdir(exist_ok=True)
    writer=Writer(folder/'independent');started=time.monotonic()
    mesh=read(out/'mesh-h0.json');ground=read(out/'ground.json');case=read(out/'stages.json')
    cfg=prior.configuration(case,g,5e-5)
    if enabled:
        cfg['joint_system']['global_step_control']={'switch_bisections':BISECTIONS}
        cfg['settings']['minimum_increment']=min(cfg['settings']['minimum_increment'],1/(g*2**(BISECTIONS+1)))
    value=reference._reference(mesh,ground,cfg,g,1,7,False,True,writer,maximum_m=5e-5,switch_bisections=BISECTIONS if enabled else None)
    writer.close();write(folder/'independent-result.json',value)
    write(folder/'independent-audit.json',dict(physical=audit(folder/'independent',mesh,ground,cfg,g),elapsed_seconds=time.monotonic()-started))
    print('independent complete',g,flush=True)


def restart(out,g,timeout):
    folder=out/f'g{g}';raw=folder/'native-raw'
    binary=out/('probe.exe' if (out/'probe.exe').exists() else 'probe')
    command=[str(binary),str(out/'mesh-h0.json'),str(out/'ground.json'),str(folder/'native-stages.json')]
    expected=without_clock(read(raw/'run.json'));selected={}
    # One real checkpoint per phase, including points after adaptive cutbacks.
    for checkpoint in sorted(raw.glob('checkpoint-*.json')):
        cp=read(checkpoint)
        state=cp['run']['last_committed']
        selected.setdefault(state['stage_index'],checkpoint)
    assert set(selected)==set(range(7)),('checkpoint phases missing',selected)
    def one(item):
        phase,checkpoint=item
        target=folder/('restart-'+checkpoint.stem);cmd=command+[str(target),'--resume',str(checkpoint)]
        started=time.monotonic()
        with (folder/(target.name+'.log')).open('x') as log: result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        assert result.returncode==0
        assert without_clock(read(target/'run.json'))==expected,'separate-process Run mismatch'
        n=int(checkpoint.stem.split('-')[-1]);count=0
        with (raw/'accepted.jsonl').open() as a,(target/'accepted.jsonl').open() as b:
            for _ in range(n): next(a)
            for x,y in itertools.zip_longest(a,b):
                assert x is not None and y is not None and json.loads(x)==json.loads(y),'restart suffix changed'
                count+=1
        row=dict(stage=STAGES[phase],checkpoint=checkpoint.name,checkpoint_sha256=sha(checkpoint),process_id=read(target/'verification.json')['process_id'],suffix_states_exact=count,entire_run_exact_except_elapsed=True,elapsed_seconds=time.monotonic()-started)
        print('restart complete',g,phase,flush=True)
        return row
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows=list(pool.map(one,selected.items()))
    target=folder/'uninterrupted';cmd=command+[str(target),'--whole']
    with (folder/'uninterrupted.log').open('x') as log: result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
    assert result.returncode==0 and without_clock(read(target/'run.json'))==expected
    write(folder/'restart-verification.json',dict(steps=g,checks=rows,uninterrupted_exact_except_elapsed=True,all_seven_phases=True))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--solver',type=Path,default=ROOT/'target/release/examples/wet_joint_step_probe')
    p.add_argument('--mode',choices=('init','native','independent','restart'),required=True);p.add_argument('--steps',type=int,nargs='+',default=STEPS)
    p.add_argument('--timeout',type=int,default=1800);p.add_argument('--legacy',action='store_true');p.add_argument('--pc-series',action='store_true')
    a=p.parse_args();out=a.output.resolve()
    if a.mode=='init': initialize(out,a.solver.resolve());return
    assert all(g in STEPS or (a.pc_series and g in (512,1024,2048)) for g in a.steps), 'long series requires explicit --pc-series on user PC'
    for g in a.steps:
        if a.mode=='native': native(out,g,a.timeout,not a.legacy)
        elif a.mode=='independent': independent(out,g,not a.legacy)
        else: restart(out,g,a.timeout)
if __name__=='__main__': main()
