"""Bounded Linux before/after experiment. A new output directory is mandatory."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'interim/2026-11-09-nc-pseudo-arclength')]
BASE = 'e278dbf17c771e87098b6180e9c65d9e4a6ad2e7'
ARC = Path('interim/2026-11-09-nc-pseudo-arclength/independent_arc.py')


def baseline(reference, output, source, summary):
    spec=importlib.util.spec_from_file_location('nc_arc_old_direct',source)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    original=np.linalg.solve; worst=None; count=0
    def observe(a,b):
        nonlocal worst,count
        x=original(a,b)
        if b[-1]==1 and np.all(b[:-1]==0):
            count+=1
            residual=float(np.linalg.norm(a@x-b)/max(np.linalg.norm(b),1e-30))
            if worst is None or residual>worst['residual']:
                worst=dict(matrix=a.tolist(),rhs=b.tolist(),solution=x.tolist(),residual=residual)
        return x
    np.linalg.solve=observe
    try:
        report=old.run(reference,output,2e-5,16 if reference.name=='reference128' else 8)
    finally:
        np.linalg.solve=original
    summary.write_text(json.dumps(dict(actual_bordered_solves=count,worst=worst),indent=2)+'\n')
    return report


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    source=output/'baseline-source.py'
    source.write_bytes(subprocess.run(['git','show',f'{BASE}:{ARC.as_posix()}'],cwd=ROOT,check=True,capture_output=True).stdout)
    commands=[]
    for steps in (128,400):
        command=[sys.executable,'tools/low_preload_nc_reference.py','--output',str(output/f'reference{steps}'),'--steps',str(steps),'--method','rk4','--budget-seconds','60']
        with (output/f'reference{steps}.log').open('w') as log:
            rc=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT).returncode
        commands.append(dict(command=command,return_code=rc))
        if rc:raise RuntimeError('reference process failed; retain all outputs')
        report=baseline(output/f'reference{steps}',output/f'before{steps}',source,output/f'before{steps}-systems-summary.json')
        (output/f'before{steps}.log').write_text(json.dumps(report,indent=2)+'\n')
    command=[sys.executable,str(ARC),'--reference',str(output/'reference400'),'--output',str(output/'after400'),'--step-size','0.00002','--count','8']
    with (output/'after400.log').open('w') as log:
        rc=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT).returncode
    commands.append(dict(command=command,return_code=rc))
    (output/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    if rc:raise RuntimeError('corrected arc failed; retain all outputs')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
