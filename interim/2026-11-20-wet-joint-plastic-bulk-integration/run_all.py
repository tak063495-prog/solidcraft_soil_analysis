"""Fresh bounded C-03 matrix; existing outputs are never overwritten."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from study import ROOT,HERE,read,write,sha,without_clock,invoke

ROUTES=[('baseline',1e-4,1e-4),('material-half',5e-5,1e-4),('material-quarter',2.5e-5,1e-4),('contact-half',1e-4,5e-5),('contact-quarter',1e-4,2.5e-5),('single-partition-control',1e-3,5e-4)]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--solver',type=Path,required=True);p.add_argument('--steps',type=int,default=32);p.add_argument('--element',choices=['tet4','tet10'],default='tet4');p.add_argument('--timeout',type=int,default=900);a=p.parse_args()
    out=a.output.resolve();solver=a.solver.resolve();assert not out.exists() and solver.is_file()
    study=HERE/'study.py';cmd=[sys.executable,str(study),'--output',str(out)]
    subprocess.run(cmd+['--mode','init','--solver',str(solver),'--cohesion','30','--open-factor','.12','--element',a.element],check=True)
    sources={}
    for directory in ['src','tools','tests','examples','interim/2026-11-20-wet-joint-plastic-bulk-integration']:
        for path in (ROOT/directory).rglob('*'):
            if path.is_file() and path.suffix in ('.py','.rs','.toml') and '__pycache__' not in path.parts:sources[str(path.relative_to(ROOT))]=sha(path)
    for path in [ROOT/'Cargo.toml',ROOT/'Cargo.lock',ROOT/'rust-toolchain.toml']:sources[path.name]=sha(path)
    write(out/'source-sha256.json',sources);processes=[]
    try:
        for name,material,contact in ROUTES:
            for mode in ('native','independent'):
                route=name+'-'+mode
                args=cmd+['--mode',mode,'--folder',route,'--steps',str(a.steps),'--material',str(material),'--contact',str(contact),'--timeout',str(a.timeout)]
                result=invoke(args,out/(route+'-driver.log'),a.timeout+60);processes.append(result)
                if result['exit_code']!=0:raise RuntimeError(('incomplete route',route,result))
        result=invoke(cmd+['--mode','restart','--folder','restarts','--source','baseline-native','--steps',str(a.steps),'--timeout',str(a.timeout)],out/'restart-driver.log',7*a.timeout+60);processes.append(result)
        if result['exit_code']!=0:raise RuntimeError('restart failed')
        result=invoke(cmd+['--mode','whole','--folder','whole-native','--steps',str(a.steps),'--material','1e-4','--contact','1e-4','--timeout',str(a.timeout)],out/'whole-driver.log',a.timeout+60);processes.append(result)
        if result['exit_code']!=0:raise RuntimeError('unobserved solve failed')
        assert without_clock(read(out/'baseline-native/raw/run.json'))==without_clock(read(out/'whole-native/raw/run.json'))
        result=invoke([sys.executable,HERE/'measure.py','--output',out,'--steps',str(a.steps),'--routes',*[x[0] for x in ROUTES],'--report',out/'measurements.json'],out/'measure-driver.log',a.timeout+60);processes.append(result)
        if result['exit_code']!=0:raise RuntimeError('full field audit failed')
        assert all(sha(ROOT/name)==digest for name,digest in sources.items()),'numerical source changed during solves'
        write(out/'completion.json',dict(complete=True,fresh_native_routes=6,fresh_independent_routes=6,accepted_states_per_route=7*a.steps,
              source_files_unchanged=len(sources),unobserved_run_exact_except_elapsed=True,real_checkpoint_restarts=7,processes=processes,
              C03_numeric_acceptance=False,claim='completed integration observations; numerical judgement is recorded separately'))
    except Exception as e:
        write(out/'matrix-failure.json',dict(complete=False,passed=False,error=repr(e),processes=processes));raise
if __name__=='__main__':main()
