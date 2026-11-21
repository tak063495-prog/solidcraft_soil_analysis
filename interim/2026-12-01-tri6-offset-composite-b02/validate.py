"""Every existing Rust/Python test, without skips or weakened checks."""
import argparse
import concurrent.futures
import os
from pathlib import Path
import sys
from cases import ROOT, write
from study import command
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[k]='1'
os.environ.setdefault('RUST_TEST_THREADS','2')
os.environ.setdefault('CARGO_BUILD_JOBS','2')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--part',choices=('rust','python'),required=True);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    def call(name,args):
        return name,command(args,out/(name+'.log'),expect=None)
    results={}
    if a.part=='rust':
        for name,cmd in [('fmt',['cargo','fmt','--all','--check']),('clippy',['cargo','clippy','--locked','--all-targets','--','-D','warnings']),
            ('build',['cargo','build','--locked','--release']),('rust-test',['cargo','test','--locked','--release','--all-targets']),
            ('rust-doc',['cargo','test','--locked','--release','--doc'])]:
            name,r=call(name,cmd);results[name]=r;write(out/'processes.json',results)
    else:
        def regular():
            rows=[]
            for group in ('python-core-0','python-core-1'):
                rows.append(call(group,[sys.executable,'tools/run_ci_python_tests.py','--group',group,'--output',out/group]))
            return rows
        def sealed():
            rows=[]
            for group in ('python-mean-evidence','python-spatial-evidence'):
                rows.append(call(group,[sys.executable,'tools/run_ci_python_tests.py','--group',group,'--output',out/group]))
            return rows
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            for rows in pool.map(lambda f:f(),(regular,sealed)):
                results.update(rows);write(out/'processes.json',results)
    return int(any(x['exit_code'] for x in results.values()))


if __name__=='__main__':sys.exit(main())
