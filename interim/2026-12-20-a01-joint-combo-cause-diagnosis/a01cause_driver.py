"""Explicit bounded input-only controls and separate native failure observations."""
import argparse
import json
from pathlib import Path
import sys
import time
import traceback
import numpy as np
from a01cause_binding import fixtures,fields,study
from a01cause_reference import Model,run


def write(path,value):
    def serial(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    with Path(path).open('x',encoding='utf-8') as f:json.dump(value,f,ensure_ascii=False,allow_nan=False,default=serial);f.write('\n')


def independent(values,out,steps,*,fixed_parts=None,boost=False,trace=False):
    out.mkdir(parents=True,exist_ok=False)
    for name,value in zip(('mesh','ground','case'),values):write(out/(name+'.json'),value)
    result=dict(complete=False,uses_native_results=False)
    started=time.monotonic()
    with (out/'fields.jsonl').open('x',encoding='utf-8') as stream,(out/'trace.jsonl').open('x',encoding='utf-8') as log:
        def observe(event):
            log.write(json.dumps(event,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v.item(),allow_nan=False)+'\n');log.flush()
        try:result=run(Model(*values),steps,fields.Sink(stream),fixed_parts=fixed_parts,boost_nontrivial=boost,observer=observe if trace else None)
        except Exception as error:result.update(error=repr(error),traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-started,uses_native_results=False,fixed_parts=fixed_parts,boost_nontrivial=boost)
    write(out/'result.json',result)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mode',choices=('controls','stacked','native'),required=True)
    p.add_argument('--steps',type=int,default=8);p.add_argument('--case',default='shared-stress_free')
    p.add_argument('--parts',type=int);p.add_argument('--boost',action='store_true')
    p.add_argument('--probe',type=Path);p.add_argument('--baseline',type=Path);a=p.parse_args()
    if a.mode=='controls':
        values=fixtures.cases(steps=a.steps)[a.case] if isinstance(fixtures.cases(steps=a.steps),dict) else None
        if values is None:raise ValueError('case inventory differs')
        result=independent(values,a.output,a.steps,fixed_parts=a.parts,boost=a.boost)
        if not result['complete']:raise RuntimeError(result['error'])
    elif a.mode=='native':
        a.output.mkdir(parents=True,exist_ok=False)
        values=fixtures.cases(steps=a.steps)[a.case]
        for name,value in zip(('mesh','ground','case'),values):write(a.output/(name+'.json'),value)
        study.native(a.output,a.probe.resolve(),a.steps,900)
    else:
        a.output.mkdir(parents=True,exist_ok=False);values=fixtures.stacked(steps=a.steps)
        for name,value in zip(('mesh','ground','case'),values):write(a.output/(name+'.json'),value)
        independent(values,a.output/'independent',a.steps,trace=True)
        for name,command in (('observed',[a.probe.resolve(),*study.inputs(a.output),a.output/'native']),
                             ('unobserved',[a.probe.resolve(),*study.inputs(a.output),a.output/'unobserved','--unobserved'])):
            process=study.invoke(command,a.output/(name+'.log'),900);write(a.output/(name+'-process.json'),process)
            if process['exit_code']!=0:raise RuntimeError(process)
        if a.baseline:
            process=study.invoke([a.baseline.resolve(),'solve-stages',*study.inputs(a.output),a.output/'baseline.json'],a.output/'baseline.log',900)
            write(a.output/'baseline-process.json',process)
        observed=fixtures.read(a.output/'native/run.json');unobserved=fixtures.read(a.output/'unobserved/run.json')
        assert study.exact(observed)==study.exact(unobserved)
        baseline=fixtures.read(a.output/'baseline.json') if a.baseline else unobserved
        import copy
        left=copy.deepcopy(observed);right=copy.deepcopy(baseline)
        # A different executable must retain its own provenance SHA. Compare
        # every numerical/history/input field, excluding only that identity.
        binary_identity={side:v['provenance'].pop('solver_binary_sha256') for side,v in (('observed',left),('baseline',right))}
        assert study.exact(left)==study.exact(right)
        write(a.output/'software-observation-equality.json',dict(exact_except_elapsed=True,baseline_binary_run_compared=bool(a.baseline),cross_binary_excluded_field='/provenance/solver_binary_sha256',binary_identity=binary_identity,numerical_complete=observed['complete'],failure_reason=observed['failure_reason']))

if __name__=='__main__':main()
