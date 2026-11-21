"""Sixteen bounded fresh global histories; P/R switches explicit, default off."""
import argparse
import itertools
import shutil
from c04follow_common import *

OLD=ROOT/'interim/2026-11-22-wet-joint-transient-water'

def specs():
    return [dict(m=m,water=32,engine=e,p=p,r=r,expected_states=7*m+192)
        for m,e,p,r in itertools.product((16,64),('native','independent'),('N','I'),('N','I'))]

def legacy():
    sys.path.insert(0,str(ROOT/'tools'));sys.path.insert(0,str(OLD))
    return load('c04follow_legacy_study',OLD/'study.py')

def worker(a):
    import c04follow_reference as reference
    s=legacy();folder=a.run/f'm{a.m}-{a.engine}-p{a.p}r{a.r}'
    if a.engine=='native':
        s.native(a.run,folder,a.m,32,3600)
        result=read(folder/'raw/run.json')
        attempts=result['geostatic']['attempts']+result['attempts']
        predictors=[q for v in attempts for q in v.get('contact_predictors',[])]
        if predictors:raise RuntimeError(('unsupported auxiliary predictor exercised',predictors))
        write(folder/'predictor-audit.json',dict(auxiliary_predictors=0,supported_case=True,attempts=len(attempts)))
    else:
        folder.mkdir();case=s.mechanical_case(a.m,32);write(folder/'case.json',case)
        model=reference.Model(read(a.run/'mesh.json'),read(a.run/'ground.json'),case)
        f,emit=sink(folder/'fields.jsonl')
        try:result=reference.run(model,a.m,case['settings']['material_maximum_substep'],case['joint_system']['maximum_increment_m'],emit,predictor=a.p,repartition=a.r)
        finally:f.close()
        write(folder/'result.json',result)
    print(folder.name, 'complete',flush=True)

def binding():
    paths=set(subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines())
    paths.update(p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.rs') or p.name=='SPEC.ja.md')
    return {p:sha(ROOT/p) for p in sorted(paths) if Path(p).suffix in ('.py','.rs','.toml','.ps1') or Path(p).name in ('Cargo.lock','SPEC.ja.md')}

def run(a):
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);build=a.build.resolve()
    write(out/'plan.json',specs());before=binding();write(out/'source-binding-before.json',before)
    for name in before:
        p=out/'source-snapshot'/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p)
    suffix='.exe' if os.name=='nt' else ''
    for name in ('mesh','ground'):shutil.copy2(OLD/('inputs/'+name+'.json'),out/(name+'.json'))
    shutil.copy2(build/('target/release/c04follow_probe'+suffix),out/('probe'+suffix))
    write(out/'identity.json',dict(base=BASE,binaries={'probe'+suffix:sha(out/('probe'+suffix))}))
    write(out/'isolated-build-binding.json',read(build/'binding.json'))
    processes=[]
    try:
        for spec in specs():
            args=[sys.executable,HERE/'c04follow_run.py','--worker','--run',out,'--m',spec['m'],'--engine',spec['engine'],'--p',spec['p'],'--r',spec['r']]
            name=f"m{spec['m']}-{spec['engine']}-p{spec['p']}r{spec['r']}"
            proc=invoke(args,out,name,{'C04FOLLOW_P':spec['p'],'C04FOLLOW_R':spec['r']})
            processes.append(dict(name=name,**proc));print(name,proc['elapsed_seconds'],flush=True)
        after=binding();write(out/'source-binding-after.json',after)
        if before!=after:raise ValueError('sources changed during numerical matrix')
        write(out/'result.json',dict(execution_complete=True,passed=False,c04_complete=False,diagnostic_only=True,processes=processes,sources_unchanged=True,paths=16,accepted_states=7552))
    except BaseException as error:
        write(out/'failure.json',dict(execution_complete=False,passed=False,error=repr(error),processes=processes));raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path);p.add_argument('--build',type=Path)
    p.add_argument('--worker',action='store_true');p.add_argument('--run',type=Path)
    p.add_argument('--engine',choices=('native','independent'));p.add_argument('--m',type=int,choices=(16,64))
    p.add_argument('--p',choices=('N','I'));p.add_argument('--r',choices=('N','I'))
    a=p.parse_args()
    if a.worker:worker(a)
    else:run(a)
