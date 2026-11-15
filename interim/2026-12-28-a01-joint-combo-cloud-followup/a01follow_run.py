"""Fresh input-only bidirectional P x R matrix and separate stacked failures."""
import argparse, copy, itertools, json, os
from pathlib import Path
import subprocess, sys, time, traceback
import numpy as np
from a01follow_binding import fixtures, fields, study, original_reference
from a01follow_policies import Predictor, serial
from a01follow_reference import Model, run

def write(path,value):
    with Path(path).open('x',encoding='utf-8') as f:
        json.dump(value,f,ensure_ascii=False,allow_nan=False,default=serial);f.write('\n')

def native(values,out,probe,predictor,repartition,steps=8):
    out.mkdir(parents=True,exist_ok=False)
    for name,v in zip(('mesh','ground','case'),values):write(out/(name+'.json'),v)
    env=dict(os.environ,A01FOLLOW_P=predictor,A01FOLLOW_R=repartition)
    started=time.monotonic()
    with (out/'process.log').open('x') as log:
        p=subprocess.Popen([str(probe),*map(str,study.inputs(out)),str(out/'native')],env=env,stdout=log,stderr=subprocess.STDOUT)
        timed_out=False
        try:p.wait(timeout=900)
        except subprocess.TimeoutExpired:timed_out=True;p.kill();p.wait()
    result=dict(exit_code=p.returncode,pid=p.pid,timed_out=timed_out,elapsed_seconds=time.monotonic()-started,predictor=predictor,repartition=repartition)
    if (out/'native/run.json').exists():
        raw=fixtures.read(out/'native/run.json');result.update(complete=raw['complete'],failure_reason=raw['failure_reason'])
        model=original_reference.Model(*values)
        names=['geostatic']+[s['name'] for s in values[2]['stages']]
        with (out/'native/accepted.jsonl').open() as a,(out/'native/constitutive.jsonl').open() as b,(out/'native/structure-reports.jsonl').open() as c,(out/'fields.jsonl').open('x') as stream:
            sink=fields.Sink(stream);count=0
            for s,d,r in itertools.zip_longest(a,b,c):
                assert s is not None and d is not None and r is not None
                state,diag,report=map(json.loads,(s,d,r));phase=count//steps;progress=(count%steps+1)/steps
                assert state['stage_index']==phase and diag['stage']==phase and diag['progress']==progress
                sink(dict(stage=names[phase],stage_index=phase,progress=progress,factor=state['boundary_factor']),fields.native_fields(model,state,diag,report));count+=1
        result['accepted_states']=count
    else:
        prefix=out/'native/accepted.jsonl'
        accepted=sum(1 for _ in prefix.open()) if prefix.exists() else 0
        result.update(complete=False,accepted_states=accepted,failure_reason='diagnostic process incomplete; partials retained')
    write(out/'result.json',result)
    return result

def independent(values,out,service,predictor='independent',repartition='independent',steps=8,boost=False):
    out.mkdir(parents=True,exist_ok=False)
    for name,v in zip(('mesh','ground','case'),values):write(out/(name+'.json'),v)
    started=time.monotonic();result=dict(complete=False);bridge=Predictor(service) if predictor=='native' else None
    with (out/'fields.jsonl').open('x') as stream,(out/'trace.jsonl').open('x') as trace:
        def observe(event):trace.write(json.dumps(event,default=serial,allow_nan=False)+'\n');trace.flush()
        try:
            result=run(Model(*values),steps,fields.Sink(stream),observer=observe,predictor=predictor,repartition=repartition,service=bridge,boost_nontrivial=boost)
        except Exception as error:result.update(error=repr(error),traceback=traceback.format_exc())
        finally:
            if bridge:bridge.close()
    result.update(pid=os.getpid(),elapsed_seconds=time.monotonic()-started,predictor=predictor,repartition=repartition,boost=boost,uses_native_fields_or_partitions=False)
    result['accepted_states']=sum(1 for _ in (out/'fields.jsonl').open())
    write(out/'result.json',result)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--case',choices=('shared-stress_free','shared-transfer_force','joint-plate','stacked'),required=True)
    p.add_argument('--engine',choices=('native','independent'),required=True);p.add_argument('--predictor',choices=('native','independent'),required=True)
    p.add_argument('--repartition',choices=('native','independent'),required=True);p.add_argument('--steps',type=int,default=8);p.add_argument('--boost',action='store_true')
    a=p.parse_args();values=fixtures.stacked(steps=a.steps) if a.case=='stacked' else fixtures.cases(steps=a.steps)[a.case]
    binary=a.build.resolve()/'target/release'
    exe=lambda name:binary/(name+('.exe' if os.name=='nt' else ''))
    if a.engine=='native':r=native(values,a.output,exe('a01follow_probe'),a.predictor,a.repartition,a.steps)
    else:r=independent(values,a.output,exe('a01follow_predictor'),a.predictor,a.repartition,a.steps,a.boost)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','traceback')},default=serial,ensure_ascii=False))

if __name__=='__main__':main()
