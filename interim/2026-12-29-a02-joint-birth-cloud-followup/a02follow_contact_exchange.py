"""P/R exchange with actual contact repartitions on the same 24 Tet10 cycle.
Uses the already studied 0.0005m contact bound, not a new physical input.
"""
import argparse,subprocess,sys,time,traceback
from pathlib import Path
from a02follow_run import worker,write
from a02follow_binding import study
MAXIMUM=.0005

def one(out,build,engine,p,r):
    original=study.fixture
    def fixture(*args,**kwargs):
        mesh,ground,case=original(*args,**kwargs)
        case['joint_system']['maximum_increment_m']=MAXIMUM
        return mesh,ground,case
    study.fixture=fixture
    try:worker(out,build,4,1,'stress_free-cycle',engine,p,r)
    finally:study.fixture=original

def batch(out,build):
    out.mkdir(parents=True,exist_ok=False);rows=[]
    for engine in ('native','independent'):
        for p in ('native','independent'):
            for r in ('native','independent'):
                name=f'{engine}-p{p}-r{r}';target=out/name;start=time.monotonic()
                with (out/(name+'.log')).open('x',encoding='utf-8') as log:
                    proc=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'worker','--output',str(target),'--build',str(build),'--engine',engine,'--predictor',p,'--repartition',r],stdout=log,stderr=subprocess.STDOUT)
                    timed=False
                    try:code=proc.wait(timeout=960)
                    except subprocess.TimeoutExpired:timed=True;proc.kill();code=proc.wait()
                row=dict(path=name,engine=engine,predictor=p,repartition=r,n=4,minimum_parts=1,maximum_contact_increment_m=MAXIMUM,pid=proc.pid,exit_code=code,timed_out=timed,elapsed_seconds=time.monotonic()-start)
                if (target/'completion.json').exists():row.update(study.read(target/'completion.json'))
                rows.append(row);write(out/(name+'.process.json'),row);print(name,code,round(row['elapsed_seconds'],3),flush=True)
    write(out/'paths.json',rows)
    if not all(r['exit_code']==0 and r.get('complete') for r in rows):raise RuntimeError('actual-contact exchange incomplete; evidence retained')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('batch','worker'));p.add_argument('--output',type=Path,required=True);p.add_argument('--build',type=Path,required=True);p.add_argument('--engine',choices=('native','independent'));p.add_argument('--predictor',choices=('native','independent'));p.add_argument('--repartition',choices=('native','independent'));a=p.parse_args()
    if a.mode=='batch':batch(a.output.resolve(),a.build.resolve())
    else:
        try:one(a.output.resolve(),a.build.resolve(),a.engine,a.predictor,a.repartition)
        except BaseException as e:
            if a.output.exists():write(a.output/'failure.json',dict(error=repr(e),traceback=traceback.format_exc(),a02_passed=False))
            raise
