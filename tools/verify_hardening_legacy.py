"""Run unchanged accepted MCC local inputs with the extended solver."""
import argparse
import json
from pathlib import Path
import hashlib
from verify_checkpoint import invoke
from verify_mcc_regression import verify


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def compare(solver,baseline,output):
    output.mkdir(parents=True,exist_ok=False)
    accepted=json.loads((baseline/'final-validation.json').read_text(encoding='utf-8-sig'))
    assert accepted['passed']
    rows=[]
    for shape in ('flat','curved'):
        folder=baseline/shape;target=output/(shape+'-result.json')
        paths=[folder/'ground/mesh.json',folder/'ground/ground.json',folder/'stages.json']
        invoke([str(solver),'solve-stages',*map(str,paths),str(target)],output/(shape+'.log'))
        report=output/(shape+'-comparison.json');verify(folder/'stage-results.json',target,report)
        rows.append({'shape':shape,'inputs_sha256':{str(p):sha(p) for p in paths},'exact_committed_states_and_frames':True})
    with (output/'verification.json').open('x') as f:json.dump({'schema':'soil-fem-hardening-legacy-regression/1','passed':True,'solver_sha256':sha(solver),'comparisons':rows},f,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();compare(a.solver.resolve(),a.baseline.resolve(),a.output_dir.resolve())
