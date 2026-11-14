"""C-01 only: bounded global-step series; measurements never certify convergence."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OLD=ROOT/'interim/2026-11-15-wet-joint-convergence'
sys.path.insert(0,str(OLD))
import study as prior
from study import read,write,sha,Writer,states,STAGES,FIELDS,fem,probe
spec=importlib.util.spec_from_file_location('c01_input_only_reference',HERE/'reference_steps.py')
probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
BASE='b069ec84da432f677c33fab5c984b6e680835bc4'
STEPS=(16,32,64,128,256)


def without_clock(v):
    if isinstance(v,dict):return {k:without_clock(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [without_clock(x) for x in v]
    return v


def initialize(out,solver):
    out.mkdir(parents=True,exist_ok=False)
    for name,source in [('mesh-h0.json','joint-water-column-tet10-mesh.json'),('ground.json','ground-joint-water-mc.json'),('stages.json','stages-joint-water-cycle.json')]:
        v=read(ROOT/'examples'/source)
        if name=='ground.json':v=fem.normalize_ground(v)
        if name=='stages.json':v=fem.normalize_stage_case(v)
        write(out/name,v)
    snapshot=out/('probe.exe' if solver.suffix=='.exe' else 'probe')
    shutil.copy2(solver,snapshot)
    write(out/'identity.json',dict(starting_commit=BASE,solver_sha256=sha(snapshot),inputs={n:sha(out/n) for n in ('mesh-h0.json','ground.json','stages.json')},physical_model_changed=False,tolerances_changed=False))


def native(out,g,timeout):
    folder=out/f'g{g}';folder.mkdir(exist_ok=True)
    mesh=read(out/'mesh-h0.json');ground=read(out/'ground.json');case=read(out/'stages.json')
    cfg=prior.configuration(case,g,5e-5)
    write(folder/'native-stages.json',cfg)
    binary=out/('probe.exe' if (out/'probe.exe').exists() else 'probe')
    cmd=[str(binary),str(out/'mesh-h0.json'),str(out/'ground.json'),str(folder/'native-stages.json'),str(folder/'native-raw')]
    started=time.monotonic()
    with (folder/'native.log').open('x') as f:
        r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
    assert r.returncode==0,('native failure',g,r.returncode)
    raw=folder/'native-raw';run=read(raw/'run.json');prior.schedule(run,g)
    assert read(raw/'verification.json')['accepted_states_written']==7*g
    op=prior.operators(mesh,ground,case,'cut');writer=Writer(folder/'native')
    with (raw/'accepted.jsonl').open() as f:
        for i,line in enumerate(f):
            s=json.loads(line);phase,step=divmod(i,g)
            writer(dict(stage=STAGES[phase],progress=(step+1)/g,factor=s['boundary_factor']),prior.native_fields(s,op))
    writer.close()
    write(folder/'native-audit.json',prior.audit(folder/'native',mesh,ground,case,'cut',g))
    write(folder/'native-timing.json',dict(elapsed_seconds=time.monotonic()-started,command=cmd,exit_code=0))
    print('native complete',g,flush=True)


def independent(out,g,mode):
    folder=out/f'g{g}';folder.mkdir(exist_ok=True)
    key={'fixed':'fixed8','cap':'independent','open-fine':'open256-rest64','open-coarse':'open64-rest256'}[mode]
    writer=Writer(folder/key);started=time.monotonic()
    mesh=read(out/'mesh-h0.json');ground=read(out/'ground.json');case=read(out/'stages.json')
    kwargs={}
    if mode=='open-fine':kwargs['progress_grids']={'open':[i/256 for i in range(1,257)]}
    if mode=='open-coarse':kwargs['progress_grids']={'open':[i/64 for i in range(1,65)]}
    value=probe._reference(mesh,ground,case,g,8 if mode=='fixed' else 1,7,False,True,writer,maximum_m=None if mode=='fixed' else 5e-5,**kwargs)
    writer.close();write(folder/(key+'-result.json'),value)
    if mode in ('fixed','cap'):write(folder/(key+'-audit.json'),prior.audit(folder/key,mesh,ground,case,'cut',g))
    write(folder/(key+'-timing.json'),dict(elapsed_seconds=time.monotonic()-started,reference_uses_native_results=False))
    print(key,'complete',g,flush=True)


def restart(out,g,timeout):
    folder=out/f'g{g}';raw=folder/'native-raw'
    binary=out/('probe.exe' if (out/'probe.exe').exists() else 'probe')
    command=[str(binary),str(out/'mesh-h0.json'),str(out/'ground.json'),str(folder/'native-stages.json')]
    expected=without_clock(read(raw/'run.json'));rows=[]
    # Each process reads an on-disk, hash-bound checkpoint. All seven phases.
    for checkpoint in sorted(raw.glob('checkpoint-*.json')):
        target=folder/('restart-'+checkpoint.stem)
        cmd=command+[str(target),'--resume',str(checkpoint)]
        started=time.monotonic()
        with (folder/(target.name+'.log')).open('x') as f:
            result=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
        assert result.returncode==0
        actual=read(target/'run.json');assert without_clock(actual)==expected,'separate-process Run mismatch'
        n=int(checkpoint.stem.split('-')[-1])
        # Also compare every subsequent committed state, not just seven endpoints.
        with (raw/'accepted.jsonl').open() as a,(target/'accepted.jsonl').open() as b:
            for _ in range(n):next(a)
            count=0
            import itertools
            for x,y in itertools.zip_longest(a,b):
                assert x is not None and y is not None and json.loads(x)==json.loads(y),'restart suffix changed'
                count+=1
        rows.append(dict(checkpoint=checkpoint.name,checkpoint_sha256=sha(checkpoint),process_id=read(target/'verification.json')['process_id'],suffix_states_exact=count,entire_run_exact_except_elapsed=True,elapsed_seconds=time.monotonic()-started))
        print('restart complete',g,n,flush=True)
    # Same executable without observation/checkpoint requests.
    target=folder/'uninterrupted';cmd=command+[str(target),'--whole']
    started=time.monotonic()
    with (folder/'uninterrupted.log').open('x') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
    assert r.returncode==0 and without_clock(read(target/'run.json'))==expected,'observer changed Run'
    write(folder/'restart-verification.json',dict(steps=g,checks=rows,uninterrupted_exact_except_elapsed=True,uninterrupted_elapsed_seconds=time.monotonic()-started,all_seven_phases=len(rows)==7))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--solver',type=Path,default=ROOT/'target/release/examples/wet_joint_step_probe')
    p.add_argument('--mode',choices=('init','native','fixed','cap','open-fine','open-coarse','restart'),required=True)
    p.add_argument('--steps',type=int,nargs='+',default=STEPS);p.add_argument('--timeout',type=int,default=1800)
    a=p.parse_args();out=a.output.resolve()
    if a.mode=='init':initialize(out,a.solver.resolve());return
    assert all(g in STEPS for g in a.steps),'longer series belongs on user PC'
    for g in a.steps:
        try:
            if a.mode=='native':native(out,g,a.timeout)
            elif a.mode=='restart':restart(out,g,a.timeout)
            else:independent(out,g,a.mode)
        except Exception as e:
            write(out/f'failure-{a.mode}-{g}.json',dict(passed=False,error=repr(e),complete=False));raise

if __name__=='__main__':main()
