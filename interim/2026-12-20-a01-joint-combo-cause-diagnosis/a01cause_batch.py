"""Finite small-case cause matrix, explicit invocation only; no CI rerun."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from a01cause_binding import ROOT,fixtures,study
from a01cause_driver import write
HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);(out/'_bin').mkdir()
    for name,path in [('probe',ROOT/'target/release/examples/joint_combo_probe'),('cause-probe',ROOT/'target/release/examples/a01_joint_cause_probe'),('linear-replay',ROOT/'target/release/examples/a01_joint_linear_replay'),('solver',ROOT/'target/release/soil-fem-3d')]:
        shutil.copy2(path,out/'_bin'/name)
    write(out/'binary-sha256.json',{p.name:study.sha(p) for p in sorted((out/'_bin').iterdir())})
    process=study.invoke([sys.executable,ROOT/'tools/freeze_sources.py',out],out/'freeze.log',120)
    if process['exit_code']!=0:raise RuntimeError(process)
    execution={p.name:study.sha(p) for p in sorted(HERE.glob('*.py'))}
    (out/'diagnostic-source').mkdir()
    for name in execution:shutil.copy2(HERE/name,out/'diagnostic-source'/name)
    write(out/'diagnostic-source/sha256.json',execution)
    ledger=[]
    def execute(name,args,timeout=900):
        result=study.invoke([sys.executable,HERE/'a01cause_driver.py','--output',out/name,*args],out/(name+'.log'),timeout)
        ledger.append(dict(name=name,**result));write(out/(name+'-process.json'),result)
    for case,steps in [('shared-stress_free',n) for n in (4,8,16)]+[('shared-transfer_force',8),('joint-plate',8)]:
        execute('native-'+case+'-g'+str(steps),['--mode','native','--case',case,'--steps',steps,'--probe',out/'_bin/probe'])
    for case in ('shared-stress_free','shared-transfer_force'):
        for steps in (4,8,16):execute('adaptive-'+case+'-g'+str(steps),['--mode','controls','--case',case,'--steps',steps])
    execute('adaptive-joint-plate-g8',['--mode','controls','--case','joint-plate','--steps',8])
    for case in ('shared-stress_free','shared-transfer_force','joint-plate'):
        execute('boost-'+case+'-g8',['--mode','controls','--case',case,'--steps',8,'--boost'])
    for parts in (1,2,4,8):execute('fixed'+str(parts)+'-shared-stress_free-g8',['--mode','controls','--case','shared-stress_free','--steps',8,'--parts',parts])
    execute('stacked',['--mode','stacked','--steps',8,'--probe',out/'_bin/cause-probe','--baseline',a.baseline.resolve()])
    result=study.invoke([sys.executable,HERE/'a01cause_linear.py','--input',out/'stacked','--output',out/'linear'],out/'linear.log',120)
    ledger.append(dict(name='linear-postprocess',**result));write(out/'linear-process.json',result)
    result=study.invoke([out/'_bin/linear-replay',out/'linear/replay-input.json',out/'linear/replay.json'],out/'linear-replay.log',120)
    ledger.append(dict(name='same-system-linear-replay',**result));write(out/'linear-replay-process.json',result)
    (out/'scope-refusals').mkdir();study.refusals(out/'scope-refusals',out/'_bin/solver',120)
    write(out/'process-ledger.json',ledger)
    # The full source manifest is frozen before the first formal solve.
    sys.path.insert(0,str(ROOT/'tools'));from freeze_sources import manifest
    frozen=fixtures.read(out/'source-snapshot.json');current=manifest()
    if frozen!=current:raise RuntimeError('sources changed during formal solve matrix')
    if execution!={p.name:study.sha(p) for p in sorted(HERE.glob('*.py'))}:raise RuntimeError('diagnostic source changed during formal matrix')
    write(out/'source-binding.json',dict(source_unchanged_during_solves=True,files=len(frozen),baseline_binary_sha256=study.sha(a.baseline)))
    print('Recorded',len(ledger),'processes; numerical A-01 acceptance remains false')

if __name__=='__main__':main()
