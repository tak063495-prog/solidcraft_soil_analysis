"""C-02 bounded input-only solves; production physics and tolerances unchanged."""
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
sys.path.insert(0,str(ROOT/'interim/2026-11-15-wet-joint-convergence'))
import study as old
from study import read,write,sha,fem,Writer,STAGES
BASE='2d2d54c48070581d088959b52984ceb348e0e85e'

def initialize(out,solver,max_level=2):
    assert solver.is_file(), 'build the observation executable first'
    out.mkdir(parents=True,exist_ok=False)
    mesh=read(ROOT/'examples/joint-water-column-tet10-mesh.json')
    for h in range(max_level+1):
        write(out/f'mesh-h{h}.json',mesh)
        if h<max_level:
            mesh,audit=old.refine(mesh);write(out/f'refinement-h{h+1}.json',audit)
    write(out/'ground.json',fem.normalize_ground(read(ROOT/'examples/ground-joint-water-mc.json')))
    write(out/'stages.json',fem.normalize_stage_case(read(ROOT/'examples/stages-joint-water-cycle.json')))
    shutil.copy2(solver,out/('probe.exe' if solver.suffix=='.exe' else 'probe'))
    binary=out/('probe.exe' if solver.suffix=='.exe' else 'probe')
    write(out/'identity.json',dict(starting_commit=BASE,solver_sha256=sha(binary),inputs={p.name:sha(p) for p in out.glob('*.json')},production_changed=False))

def native(out,h,g,timeout):
    folder=out/f'h{h}-g{g}';folder.mkdir(exist_ok=True)
    binary=out/('probe.exe' if (out/'probe.exe').exists() else 'probe')
    assert sha(binary)==read(out/'identity.json')['solver_sha256']
    cfg=old.configuration(read(out/'stages.json'),g,5e-5)
    write(folder/'native-stages.json',cfg)
    cmd=[str(binary),str(out/f'mesh-h{h}.json'),str(out/'ground.json'),str(folder/'native-stages.json'),str(folder/'native-raw')]
    started=time.monotonic()
    with (folder/'native.log').open('x') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
    assert r.returncode==0,('incomplete',h,g,r.returncode)
    convert_native(out,h,g)
    schedule=old.schedule(read(folder/'native-raw/run.json'),g)
    write(folder/'native-timing.json',dict(elapsed_seconds=time.monotonic()-started,exit_code=0,command=cmd,schedule=schedule))
    print('native complete',h,g,flush=True)

def convert_native(out,h,g):
    folder=out/f'h{h}-g{g}'
    raw=folder/'native-raw';run=read(raw/'run.json');schedule=old.schedule(run,g)
    assert read(raw/'verification.json')['accepted_states_written']==7*g
    mesh=read(out/f'mesh-h{h}.json');ground=read(out/'ground.json');case=read(out/'stages.json')
    op=old.operators(mesh,ground,case,'cut');writer=Writer(folder/'native')
    enums=[]
    with (raw/'accepted.jsonl').open() as stream:
        for i,line in enumerate(stream):
            state=json.loads(line);phase,step=divmod(i,g)
            writer(dict(stage=STAGES[phase],progress=(step+1)/g,factor=state['boundary_factor']),old.native_fields(state,op))
            enums.append([p['status'] for p in state['joints']['points']])
    writer.close();write(folder/'endpoint-states.json',enums)
    write(folder/'native-audit.json',old.audit(folder/'native',mesh,ground,case,'cut',g))

def independent(out,h,g,timeout):
    folder=out/f'h{h}-g{g}';folder.mkdir(exist_ok=True)
    started=time.monotonic();writer=Writer(folder/'independent')
    mesh=read(out/f'mesh-h{h}.json');ground=read(out/'ground.json');case=read(out/'stages.json')
    result=old.probe._reference(mesh,ground,case,g,1,7,False,True,writer,maximum_m=5e-5,maximum_dofs=15000)
    writer.close();write(folder/'independent-result.json',result)
    write(folder/'independent-audit.json',old.audit(folder/'independent',mesh,ground,case,'cut',g))
    write(folder/'independent-timing.json',dict(elapsed_seconds=time.monotonic()-started,reference_uses_native_fields=False))
    print('independent complete',h,g,flush=True)

def without_clock(v):
    if isinstance(v,dict):return {k:without_clock(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [without_clock(x) for x in v]
    return v

def restart(out,h,g,timeout):
    folder=out/f'h{h}-g{g}';raw=folder/'native-raw'
    binary=out/('probe.exe' if (out/'probe.exe').exists() else 'probe')
    assert sha(binary)==read(out/'identity.json')['solver_sha256']
    expected=without_clock(read(raw/'run.json'));checks=[]
    # Use last phase midpoint at the finest spatial level: actual on-disk checkpoint.
    cp=sorted(raw.glob('checkpoint-*.json'))[-1]
    target=folder/('restart-'+cp.stem)
    cmd=[str(binary),str(out/f'mesh-h{h}.json'),str(out/'ground.json'),str(folder/'native-stages.json'),str(target),'--resume',str(cp)]
    started=time.monotonic()
    with (folder/'restart.log').open('x') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
    assert r.returncode==0 and without_clock(read(target/'run.json'))==expected
    import itertools
    with (raw/'accepted.jsonl').open() as a,(target/'accepted.jsonl').open() as b:
        for _ in range(int(cp.stem.split('-')[-1])):next(a)
        count=0
        for x,y in itertools.zip_longest(a,b):
            assert x is not None and y is not None and json.loads(x)==json.loads(y)
            count+=1
    checks.append(dict(checkpoint=cp.name,sha256=sha(cp),suffix_states_exact=count,entire_run_exact_except_elapsed=True,process_id=read(target/'verification.json')['process_id']))
    write(folder/'restart-verification.json',dict(checks=checks,finest_spatial_level=h,steps=g,solver_sha256=sha(binary),elapsed_seconds=time.monotonic()-started,real_file_separate_process=True))

def continue_native(out,h,g,timeout):
    """A timed-out attempt remains failed. Complete its path from a real CP."""
    folder=out/f'h{h}-g{g}';raw=folder/'native-raw'
    failures=list(out.glob(f'failure-native-h{h}-g{g}-*.json'))
    assert failures and all(not read(p)['complete'] for p in failures)
    assert not (raw/'run.json').exists(), 'do not overwrite a completed solve'
    saved=folder/'native-interrupted-raw';raw.rename(saved)
    cp=sorted(saved.glob('checkpoint-*.json'))[-1];n=int(cp.stem.split('-')[-1])
    binary=out/('probe.exe' if (out/'probe.exe').exists() else 'probe')
    assert sha(binary)==read(out/'identity.json')['solver_sha256']
    suffix=folder/'native-resumed-raw';started=time.monotonic()
    cmd=[str(binary),str(out/f'mesh-h{h}.json'),str(out/'ground.json'),str(folder/'native-stages.json'),str(suffix),'--resume',str(cp)]
    with (folder/'continuation.log').open('x') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
    assert r.returncode==0
    raw.mkdir();overlap=0
    with (saved/'accepted.jsonl').open() as a,(suffix/'accepted.jsonl').open() as b,(raw/'accepted.jsonl').open('x') as result:
        for _ in range(n):result.write(next(a))
        old_lines=list(a)
        for i,line in enumerate(b):
            if i<len(old_lines):assert json.loads(line)==json.loads(old_lines[i]);overlap+=1
            result.write(line)
    for name in ['run.json']:
        shutil.copy2(suffix/name,raw/name)
    for p in saved.glob('checkpoint-*.json'):shutil.copy2(p,raw/p.name)
    write(raw/'verification.json',dict(complete=True,accepted_states_written=7*g,process_id=read(suffix/'verification.json')['process_id'],assembled_from_explicit_checkpoint_continuation=True))
    write(raw/'continuation-provenance.json',dict(original_attempt_failed=True,original_failure_files=[p.name for p in failures],checkpoint=str(cp.relative_to(out)),checkpoint_sha256=sha(cp),prefix_states=n,overlap_states_exact=overlap,suffix_states=read(suffix/'verification.json')['accepted_states_written'],separate_process=True,solver_sha256=sha(binary),elapsed_seconds=time.monotonic()-started))
    convert_native(out,h,g)
    write(folder/'native-timing.json',dict(elapsed_seconds=time.monotonic()-started,exit_code=0,command=cmd,continuation_only=True,original_timeout_not_success=True,schedule=old.schedule(read(raw/'run.json'),g)))
    print('checkpoint continuation complete',h,g,flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=['init','native','independent','restart','continue-native'],required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--solver',type=Path,default=ROOT/'target/release/examples/wet_joint_spatial_probe')
    p.add_argument('--h',type=int,default=0);p.add_argument('--max-level',type=int,default=2);p.add_argument('--user-pc',action='store_true');p.add_argument('--steps',type=int,default=32);p.add_argument('--timeout',type=int,default=1800)
    a=p.parse_args();out=a.output.resolve()
    assert a.h>=0 and a.max_level>=2 and a.steps>=2
    assert a.user_pc or max(a.h,a.max_level)<=2, 'finer spatial levels belong on the user PC'
    if a.user_pc:
        import os
        os.environ['C02_USER_PC']='1'
    try:
        if a.mode=='init':initialize(out,a.solver.resolve(),a.max_level)
        else:globals()[a.mode.replace('-','_')](out,a.h,a.steps,a.timeout)
    except Exception as e:
        if out.exists():write(out/f'failure-{a.mode}-h{a.h}-g{a.steps}-{time.time_ns()}.json',dict(passed=False,error=repr(e),complete=False))
        raise
if __name__=='__main__':main()
