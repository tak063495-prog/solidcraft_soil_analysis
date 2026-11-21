"""Run the declared focused checks; record failures and actual test counts."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent

def checks(out,evidence):
    out.mkdir(parents=True,exist_ok=False);rows=[]
    commands=[('root-fmt',['cargo','fmt','--all','--check']),
        ('root-linalg',['cargo','test','--locked','--release','--lib','linalg::tests']),
        ('root-midpoint',['cargo','test','--locked','--release','--test','mcc_midpoint']),
        ('root-free-contact',['cargo','test','--locked','--release','--test','free_contact']),
        ('isolated-linear',['cargo','test','--locked','--release','--manifest-path',evidence/'experiment/retain/Cargo.toml',
            '--target-dir',evidence/'experiment/target','--test','nc_linear_experiment']),
        ('isolated-clippy',['cargo','clippy','--locked','--release','--manifest-path',evidence/'experiment/retain/Cargo.toml',
            '--target-dir',evidence/'experiment/target','--example','linear_probe','--example','nc_stop_probe',
            '--test','nc_linear_experiment','--','-D','warnings']),
        ('new-audit-guards',[sys.executable,HERE/'test_audit.py','--evidence',evidence,'-v'])]
    for name in ('test_nc_midpoint_resume.py','test_low_preload_nc_diagnosis.py','test_nc_newton_trace.py','test_low_preload_nc_followup.py'):
        commands.append((name,[sys.executable,'-m','unittest','discover','-s','tests','-p',name,'-v']))
    for label,args in commands:
        args=list(map(str,args));start=time.monotonic();path=out/(label+'.log')
        with path.open('x') as log:
            r=subprocess.run(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=240)
        content=path.read_text();rust=re.findall(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',content)
        python=re.search(r'Ran (\d+) tests?',content)
        row=dict(name=label,command=args,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start,
            rust_summaries=rust,python_tests=None if python is None else int(python.group(1)))
        rows.append(row);(out/'processes.json').write_text(json.dumps(rows,indent=2)+'\n')
    return all(r['exit_code']==0 for r in rows)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--evidence',type=Path,required=True)
    a=p.parse_args();sys.exit(0 if checks(a.output.resolve(),a.evidence.resolve()) else 1)
