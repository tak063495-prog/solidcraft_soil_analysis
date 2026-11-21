"""Explicit small matrix. Raw files are never overwritten after failures."""
import argparse
import itertools
import shutil
from a04follow_common import *

def specs():
    return dict(base=BASE,diagnostic_only=True,planned_before_solves=True,a04_complete=False,
                soil_cross=[dict(layers=h,n=n,order=q,engine=eng) for h,n,q in itertools.product((2,4,8),(4,8,16),(3,5,7)) for eng in (('native','independent') if q==3 else ('independent',))],
                order=[dict(kind=k,e=e,f=f,q=q,engine=eng) for k in ('soil','dry','hydro') for e,f,q in itertools.product((False,True),repeat=3) for eng in ('native','independent')],
                exchange=[dict(engine=eng,p=p,r=r) for eng,p,r in itertools.product(('native','independent'),('N','I'),('N','I'))])

def snapshot(out):
    scopes=('src','tools','examples','interim/2026-11-20-wet-joint-plastic-bulk-integration','interim/2026-11-22-wet-joint-transient-water','interim/2026-11-25-joint-birth-projection-a02','interim/2026-11-26-wet-soil-birth-darcy-a03','interim/2026-11-29-richards-wet-joint-a04')
    paths=[p.relative_to(ROOT).as_posix() for scope in scopes for p in (ROOT/scope).rglob('*') if p.is_file() and p.suffix in ('.rs','.py','.toml','.ps1') and '__pycache__' not in p.parts]
    paths += ['Cargo.toml','Cargo.lock','rust-toolchain.toml']
    paths+= [p.relative_to(ROOT).as_posix() for p in HERE.glob('*') if p.suffix=='.py' or p.name=='SPEC.ja.md']
    binding={p:sha(ROOT/p) for p in sorted(set(paths))}
    write(out/'source-binding.json',binding)
    for p in binding:
        target=out/'source-snapshot'/p;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/p,target)
    return binding

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--part',choices=('soil','order','exchange','all'),default='all');a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);build=a.build.resolve()
    plan=specs();write(out/'plan.json',plan);binding=snapshot(out);write(out/'isolated-build-binding.json',read(build/'binding.json'))
    suffix='.exe' if os.name=='nt' else '';binpath=build/'target/release'
    solver=binpath/('soil-fem-3d'+suffix);processes=[]
    def run(args,folder,name,switch=None):
        result=invoke(args,folder,name,switch);processes.append(dict(name=str(folder.relative_to(out))+'/'+name,**result));print(processes[-1]['name'],result['elapsed_seconds'],flush=True)
    try:
        if a.part in ('soil','all'):
            from a04follow_soil import fixture
            for spec in plan['soil_cross']:
                h,n,q,eng=[spec[k] for k in ('layers','n','order','engine')]
                folder=out/'soil'/f'h{h}-n{n}-q{q}'
                if not folder.exists():
                    folder.mkdir(parents=True)
                    for k,v in zip(('mesh','ground','case'),fixture(h,n)):write(folder/(k+'.json'),v)
                if eng=='native':
                    run([solver,'solve-one-way-construction',*[folder/(k+'.json') for k in ('mesh','ground','case')],folder/'native.json'],folder,'native-solve')
                run([sys.executable,HERE/'a04follow_soil.py','--folder',folder,'--engine',eng,'--order',q],folder,eng+'-fields')
        if a.part in ('order','all'):
            from a04follow_soil import fixture
            for spec in plan['order']:
                kind,e,f,q,eng=[spec[k] for k in ('kind','e','f','q','engine')]
                folder=out/'order'/kind/f'e{int(e)}f{int(f)}q{int(q)}'
                switch={'A04FOLLOW_Q':'1'} if q else {}
                if kind=='soil':
                    if not folder.exists():
                        folder.mkdir(parents=True)
                        for k,v in zip(('mesh','ground','case'),fixture(2,4,e,f)):write(folder/(k+'.json'),v)
                    if eng=='native':run([solver,'solve-one-way-construction',*[folder/(k+'.json') for k in ('mesh','ground','case')],folder/'native.json'],folder,'native-solve',switch)
                    run([sys.executable,HERE/'a04follow_soil.py','--folder',folder,'--engine',eng]+(['--reverse'] if q else []),folder,eng+'-fields',switch)
                else:
                    # Worker creates the folder, so process logs live in its parent.
                    parent=folder.parent;parent.mkdir(parents=True,exist_ok=True)
                    args=[sys.executable,HERE/'a04follow_joint_order.py','--folder',folder,'--engine',eng,'--kind',kind,'--probe',binpath/('a04follow_'+kind+'_probe'+suffix),'--solver',solver]
                    args+=sum((['--'+k] if v else [] for k,v in (('e',e),('f',f),('q',q))),[])
                    run(args,parent,folder.name+'-'+eng,switch)
        if a.part in ('exchange','all'):
            folder=out/'exchange';folder.mkdir()
            old=ROOT/'interim/2026-11-22-wet-joint-transient-water'
            for k in ('mesh','ground'):shutil.copy2(old/('inputs/'+k+'.json'),folder/(k+'.json'))
            for name,source in (('probe','a04follow_c04_probe'),('solver','soil-fem-3d')):shutil.copy2(binpath/(source+suffix),folder/(name+suffix))
            write(folder/'identity.json',dict(binaries={name+suffix:sha(folder/(name+suffix)) for name in ('probe','solver')}))
            for spec in plan['exchange']:
                eng,pp,rr=[spec[k] for k in ('engine','p','r')]
                run([sys.executable,HERE/'a04follow_c04.py','--run',folder,'--engine',eng,'--p',pp,'--r',rr],folder,eng+'-p'+pp+'r'+rr,{'A04FOLLOW_P':pp,'A04FOLLOW_R':rr})
        changed=[p for p,d in binding.items() if sha(ROOT/p)!=d]
        if changed:raise RuntimeError(('sources changed during solves',changed))
        write(out/'result.json',dict(execution_complete=True,passed=False,a04_complete=False,diagnostic_only=True,processes=processes,sources_unchanged=True,part=a.part))
    except BaseException as error:
        write(out/'failure.json',dict(execution_complete=False,passed=False,error=repr(error),processes=processes));raise

if __name__=='__main__':main()
