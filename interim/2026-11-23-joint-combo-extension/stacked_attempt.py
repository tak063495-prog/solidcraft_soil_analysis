"""Retain an interacting two-plane attempt, including numerical noncompletion."""
import argparse
from pathlib import Path
import traceback
from fixtures import ROOT,stacked,write,read
from reference import Model,run
from fields import Sink
from study import inputs,invoke


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=8)
    p.add_argument('--solver',type=Path,default=ROOT/'target/release/soil-fem-3d');p.add_argument('--timeout',type=int,default=900);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    for path,value in zip(inputs(out),stacked(steps=a.steps)):write(path,value)
    independent=dict(complete=False)
    try:
        with (out/'independent-fields.jsonl').open('x') as stream:independent=run(Model(*(read(x) for x in inputs(out))),a.steps,Sink(stream))
    except Exception as error:independent.update(error=repr(error),traceback=traceback.format_exc())
    write(out/'independent-result.json',independent)
    process=invoke([a.solver.resolve(),'solve-stages',*inputs(out),out/'native.json'],out/'native.log',a.timeout)
    write(out/'native-process.json',process)
    native=read(out/'native.json') if (out/'native.json').exists() else {}
    write(out/'verification.json',dict(passed=False,A01='未達',case='two parallel disjoint interfaces sharing a free middle soil block',
        native_complete=native.get('complete',False),native_failure_reason=native.get('failure_reason'),
        independent_complete=independent['complete'],independent_failure=independent.get('error'),process=process,
        interpretation='Numerical noncompletion is retained, not a schema refusal, proof of physical failure, or successful full-history comparison. No added supports, relaxed tolerance, or solver-default change.'))


if __name__=='__main__':main()
