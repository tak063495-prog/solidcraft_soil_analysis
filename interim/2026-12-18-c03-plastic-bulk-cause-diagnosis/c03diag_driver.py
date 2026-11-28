"""Explicit opt-in, bounded C-03 experiments. Fresh folders; failures retained."""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from c03diag_study import (HERE, ROOT, initialize, native, independent, restart,
                          invoke, configuration, Model, read, write, sha, without_clock)
from c03diag_reference import run


def routes(steps):
    result = []
    for mi, material in enumerate((1e-4, 5e-5, 2.5e-5)):
        for ci, contact in enumerate((1e-4, 5e-5, 2.5e-5)):
            result.append(dict(name=f'm{mi}c{ci}', steps=steps, material=material, contact=contact))
    result += [dict(name='m0c3', steps=steps, material=1e-4, contact=1.25e-5),
               dict(name='m2c3', steps=steps, material=2.5e-5, contact=1.25e-5),
               dict(name='g2-m2c2', steps=2*steps, material=2.5e-5, contact=2.5e-5),
               dict(name='single', steps=steps, material=1e-3, contact=5e-4)]
    return result


def fixed_reference(out, spec):
    folder = out/spec['name']; folder.mkdir()
    case = configuration(out, spec['steps'], spec['material'], spec['contact'])
    write(folder/'stages.json', case)
    model = Model(read(out/'mesh.json'), read(out/'ground.json'), case)
    started = time.monotonic()
    with (folder/'fields.jsonl').open('x', encoding='utf-8') as stream:
        def sink(row, fields):
            stream.write(json.dumps(dict(row=row, fields={k:v.tolist() for k,v in fields.items()}), allow_nan=False)+'\n')
            stream.flush()
        result = run(model, spec['steps'], spec['material'], spec['contact'], sink,
                     fixed_material=spec['fixed_material'], fixed_contact=spec['fixed_contact'])
    write(folder/'result.json', result)
    write(folder/'process.json', dict(pid=os.getpid(), elapsed_seconds=time.monotonic()-started,
                                      reference_uses_native_results=False, fixed_material=spec['fixed_material'],
                                      fixed_contact=spec['fixed_contact']))


def matrix(args):
    out=args.output.resolve()
    initialize(out,args.solver.resolve(),args.element,30.,open_factor=.12)
    specs=[]
    for spec in routes(args.steps):
        for mode in ('native','independent'):
            specs.append(dict(spec, name=spec['name']+'-'+mode, mode=mode))
    for m in (8,32):
        for c in (4,16):
            specs.append(dict(name=f'fixed-m{m}c{c}-independent',mode='fixed',steps=args.steps,
                              material=1e-4,contact=1e-4,fixed_material=m,fixed_contact=c))
    write(out/'specifications.json',specs)
    processes=[]
    for spec in specs:
        command=[sys.executable,HERE/'c03diag_driver.py','--mode','route','--output',out,
                 '--spec',json.dumps(spec),'--timeout',str(args.timeout)]
        process=invoke(command,out/(spec['name']+'-driver.log'),args.timeout+30)
        processes.append(dict(spec=spec,process=process))
        write(out/(spec['name']+'-terminal.json'),process)
    # Same executable, true files and distinct PIDs; whole run observes nothing.
    for mode in ('restart','whole'):
        command=[sys.executable,HERE/'c03diag_driver.py','--mode',mode,'--output',out,
                 '--steps',str(args.steps),'--timeout',str(args.timeout)]
        processes.append(dict(mode=mode,process=invoke(command,out/(mode+'-driver.log'),7*args.timeout+30)))
    write(out/'matrix-terminal.json',dict(processes=processes,passed=False,
          complete=all(p['process']['exit_code']==0 for p in processes),pc_executed=False))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['matrix','route','restart','whole'],required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--solver',type=Path)
    p.add_argument('--element',choices=['tet4','tet10'],default='tet4')
    p.add_argument('--steps',type=int,default=32);p.add_argument('--timeout',type=int,default=900)
    p.add_argument('--spec');args=p.parse_args();out=args.output.resolve()
    try:
        if args.mode=='matrix':matrix(args)
        elif args.mode=='restart':restart(out,out/'m0c0-native',out/'restarts',args.timeout,args.steps)
        elif args.mode=='whole':
            native(out,out/'whole',args.steps,1e-4,1e-4,args.timeout,True)
            assert without_clock(read(out/'m0c0-native/raw/run.json'))==without_clock(read(out/'whole/raw/run.json'))
        else:
            spec=json.loads(args.spec)
            if spec['mode']=='fixed':fixed_reference(out,spec)
            elif spec['mode']=='native':native(out,out/spec['name'],spec['steps'],spec['material'],spec['contact'],args.timeout)
            else:independent(out,out/spec['name'],spec['steps'],spec['material'],spec['contact'])
    except Exception as exc:
        if out.exists():write(out/f'failure-{args.mode}-{time.time_ns()}.json',dict(complete=False,passed=False,error=repr(exc)))
        raise

if __name__=='__main__':main()
