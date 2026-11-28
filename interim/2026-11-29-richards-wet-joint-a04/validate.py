"""Run all existing Rust/Python tests without skips or changed tolerances."""
import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[2]
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
os.environ.setdefault('RUST_TEST_THREADS','2');os.environ.setdefault('CARGO_BUILD_JOBS','2')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    def call(name,command):
        started=time.monotonic()
        with (out/(name+'.log')).open('x') as f:
            f.write(repr(command)+'\n');f.flush();child=subprocess.Popen(command,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT);code=child.wait()
        record=dict(command=command,pid=child.pid,exit_code=code,elapsed_seconds=time.monotonic()-started)
        (out/(name+'-process.json')).write_text(json.dumps(record,indent=2)+'\n');return record
    def rust():
        results={}
        for name,command in [('fmt',['cargo','fmt','--all','--check']),('clippy',['cargo','clippy','--locked','--all-targets','--','-D','warnings']),('rust-test',['cargo','test','--locked','--release']),('rust-examples',['cargo','test','--locked','--release','--examples'])]:
            # Continue independent checks even if an earlier one failed.
            results[name]=call(name,command)
        return results
    def cores():
        results={}
        for group in ('python-core-0','python-core-1'):
            results[group]=call(group,[sys.executable,'tools/run_ci_python_tests.py','--group',group,'--output',str(out/group)])
        return results
    def sealed():
        results={}
        for group in ('python-mean-evidence','python-spatial-evidence'):
            results[group]=call(group,[sys.executable,'tools/run_ci_python_tests.py','--group',group,'--output',str(out/group)])
        return results
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(f) for f in (rust,cores,sealed)];results={}
        for future in concurrent.futures.as_completed(futures):results.update(future.result());(out/'processes.json').write_text(json.dumps(results,indent=2)+'\n')
    (out/'processes.json').write_text(json.dumps(results,indent=2)+'\n')
    return int(any(r['exit_code'] for r in results.values()))
if __name__=='__main__':sys.exit(main())
