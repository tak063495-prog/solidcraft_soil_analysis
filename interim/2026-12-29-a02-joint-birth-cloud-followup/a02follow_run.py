"""Explicit opt-in, 24 Tet10 only; fresh solves with source/binary bindings."""
import argparse,hashlib,json,os,subprocess,sys,time,traceback
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from a02follow_binding import ROOT,HERE,study,fields,reference
from a02follow_policies import Predictor,serial
import a02follow_reference as controlled
CASES=[('transfer_force','transfer_force',False),('stress_free','stress_free',False),('stress_free-cycle','stress_free',True)]
def write(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,ensure_ascii=False,indent=2,allow_nan=False,default=serial);f.write('\n')
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def worker(out,build,n,parts,case,engine,predictor,repartition):
    name,policy,cycle=next(c for c in CASES if c[0]==case)
    values=study.fixture(policy,n,cycle);assert len(values[0]['tetrahedra'])==24 and len(values[0]['nodes'])==81
    study.make_inputs(out,values)
    binary=build/'target/release'/('a02follow_probe'+('.exe' if os.name=='nt' else ''))
    binding=dict(engine=engine,predictor=predictor,repartition=repartition,minimum_parts=parts,steps=n,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(HERE.glob('*.py'))},
        build_binding=study.read(build/'binding.json'),input_sha256={p.name:sha(p) for p in study.input_paths(out)})
    write(out/'binding.json',binding)
    if engine=='native':
        os.environ.update(A02FOLLOW_P=predictor,A02FOLLOW_R=repartition,A02FOLLOW_MINIMUM=str(parts),A02FOLLOW_TRACE="1")
        result=study.native(out,binary,n)
        states=len(list(fields.load_fields(out/'native-fields.jsonl')))
    else:
        service=Predictor(build/'target/release'/('a02follow_predictor'+('.exe' if os.name=='nt' else ''))) if predictor=='native' else None
        try:
            with (out/'independent-fields.jsonl').open('x',encoding='utf-8',newline='\n') as f,(out/'events.jsonl').open('x',encoding='utf-8',newline='\n') as events:
                def observe(v):
                    # Matrix service is input-only; preserve compact branch/repartition observations.
                    v={k:x for k,x in v.items() if k in ('kind','stage','progress','iteration','relative','counts','gap','traction','target','reset','cut','used','trust','accepted','projected_load','original_linear_relative')}
                    events.write(json.dumps(v,default=serial,allow_nan=False)+'\n');events.flush()
                result=controlled.run(controlled.Model(*values),n,fields.Sink(f),minimum_parts=parts,
                    predictor=predictor,repartition=repartition,service=service,observer=observe)
            write(out/'independent-result.json',result);states=result['accepted_states']
        finally:
            if service:service.close()
    write(out/'completion.json',dict(complete=result['complete'],accepted_states=states,a02_passed=False))
def batch(out,build):
    out.mkdir(parents=True,exist_ok=False);jobs=[]
    for n in (2,4,8):
        for parts in (2,4,8):
            for case,_,_ in CASES:
                for engine in ('native','independent'):
                    jobs.append((f'cross/n{n}-m{parts}/{case}/{engine}',n,parts,case,engine,engine,engine))
    for case,_,_ in CASES:
        for engine in ('native','independent'):
            for p in ('native','independent'):
                for r in ('native','independent'):
                    jobs.append((f'exchange/{case}/{engine}-p{p}-r{r}',4,1,case,engine,p,r))
    results=[]
    for path,n,parts,case,engine,p,r in jobs:
        target=out/path;target.parent.mkdir(parents=True,exist_ok=True)
        command=[sys.executable,str(Path(__file__).resolve()),'worker','--output',str(target),'--build',str(build),'--n',str(n),'--parts',str(parts),'--case',case,'--engine',engine,'--predictor',p,'--repartition',r]
        start=time.monotonic()
        with target.with_suffix('.log').open('x',encoding='utf-8') as log:
            proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
            timed_out=False
            try:code=proc.wait(timeout=960)
            except subprocess.TimeoutExpired:timed_out=True;proc.kill();code=proc.wait()
        row=dict(path=path,n=n,minimum_parts=parts,case=case,engine=engine,predictor=p,repartition=r,pid=proc.pid,exit_code=code,timed_out=timed_out,elapsed_seconds=time.monotonic()-start)
        if (target/'completion.json').exists():row.update(study.read(target/'completion.json'))
        results.append(row);print(path,code,round(row['elapsed_seconds'],3),flush=True)
        write(target.with_suffix('.process.json'),row)
    write(out/'paths.json',results)
    if any(r['exit_code'] or not r.get('complete') for r in results):raise RuntimeError('Incomplete diagnostic paths retained; never scope refusal or physical failure')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('batch','worker'));p.add_argument('--output',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--n',type=int,choices=(2,4,8));p.add_argument('--parts',type=int,choices=(1,2,4,8));p.add_argument('--case',choices=[x[0] for x in CASES]);p.add_argument('--engine',choices=('native','independent'));p.add_argument('--predictor',choices=('native','independent'));p.add_argument('--repartition',choices=('native','independent'));a=p.parse_args()
    if a.mode=='batch':batch(a.output.resolve(),a.build.resolve())
    else:
        try:worker(a.output.resolve(),a.build.resolve(),a.n,a.parts,a.case,a.engine,a.predictor,a.repartition)
        except BaseException as e:
            if a.output.exists():write(a.output/'failure.json',dict(error=repr(e),traceback=traceback.format_exc(),a02_passed=False))
            raise
