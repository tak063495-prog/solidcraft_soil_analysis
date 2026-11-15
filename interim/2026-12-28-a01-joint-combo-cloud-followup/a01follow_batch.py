"""Bounded fresh experiment matrix; failed paths retain all files."""
import argparse, concurrent.futures, hashlib, itertools, json, os
from pathlib import Path
import subprocess, sys, time
from a01follow_binding import ROOT
from a01follow_run import write

HERE=Path(__file__).resolve().parent

def binding(build):
    files=[HERE/name for name in ('SPEC.ja.md','a01follow_binding.py','a01follow_build.py','a01follow_service.rs',
        'a01follow_policies.py','a01follow_reference.py','a01follow_run.py','a01follow_batch.py')]
    files += sorted((ROOT/'src').rglob('*.rs'))
    files += [build/'target/release'/(name+('.exe' if os.name=='nt' else '')) for name in ('a01follow_probe','a01follow_predictor')]
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--stacked-offdiagonal',action='store_true');a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);build=a.build.resolve()
    write(out/'source-binding-before.json',binding(build))
    tasks=[]
    for case,engine,pred,part in itertools.product(('shared-stress_free','shared-transfer_force','joint-plate'),('native','independent'),('native','independent'),('native','independent')):
        tasks.append((case,engine,pred,part,8,False))
    tasks += [(case,'independent','independent','independent',8,True) for case in ('shared-stress_free','shared-transfer_force','joint-plate')]
    tasks += [('stacked',engine,pred,pred,8,False) for engine in ('native','independent') for pred in ('native','independent')]
    tasks += [('stacked','independent','independent','independent',n,False) for n in (4,16)]
    if a.stacked_offdiagonal:tasks=[('stacked',engine,pred,part,8,False) for engine,pred,part in itertools.product(('native','independent'),repeat=3) if pred!=part]
    def execute(task):
        case,engine,pred,part,steps,boost=task
        name=f'{case}__{engine}__P{pred}__R{part}__g{steps}'+('__boost' if boost else '')
        args=[sys.executable,str(HERE/'a01follow_run.py'),'--output',str(out/name),'--build',str(build),
            '--case',case,'--engine',engine,'--predictor',pred,'--repartition',part,'--steps',str(steps)]
        if boost:args+=['--boost']
        started=time.monotonic()
        with (out/(name+'.log')).open('x') as log:
            p=subprocess.Popen(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            timed_out=False
            try:code=p.wait(timeout=900)
            except subprocess.TimeoutExpired:
                timed_out=True;p.kill();code=p.wait()
                folder=out/name
                if not (folder/'result.json').exists():
                    states=sum(1 for _ in (folder/'fields.jsonl').open()) if (folder/'fields.jsonl').exists() else 0
                    write(folder/'result.json',dict(complete=False,diagnostic_only=True,pid=p.pid,
                        predictor=pred,repartition=part,accepted_states=states,elapsed_seconds=time.monotonic()-started,
                        failure_reason='900-second diagnostic wall-clock budget exhausted; partials retained',
                        numerical_outcome='incomplete',scope_refusal=False,physical_failure=False))
        result=dict(name=name,task=task,pid=p.pid,command=args,exit_code=code,timed_out=timed_out,elapsed_seconds=time.monotonic()-started)
        write(out/(name+'-process.json'),result)
        print(name,code,flush=True)
        return result
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(execute,tasks))
    write(out/'processes.json',results);write(out/'source-binding-after.json',binding(build))
    assert all(r['exit_code']==0 or r['timed_out'] for r in results)
    assert binding(build)==json.loads((out/'source-binding-before.json').read_text())

if __name__=='__main__':main()
