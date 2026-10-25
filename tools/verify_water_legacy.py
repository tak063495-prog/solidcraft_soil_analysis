"""Exact unchanged MCC local and IDH column states after adding staged water."""
import argparse
import hashlib
import json
from pathlib import Path
from verify_checkpoint import invoke
from verify_mcc_regression import verify


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def compare(solver,output):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False);rows=[]
    baselines={'mcc-box':root/'runs/local-acceptance-20261024-033221-8ada5d9c/flat',
        'mcc-cylinder':root/'runs/local-acceptance-20261024-033221-8ada5d9c/curved',
        'idh-tet4':root/'runs/hardening-acceptance-20261024-073348-c936ff06/tet4',
        'idh-tet10':root/'runs/hardening-acceptance-20261024-073348-c936ff06/tet10'}
    for name,folder in baselines.items():
        paths=[folder/'ground/mesh.json',folder/'ground/ground.json',folder/'stages.json'];target=output/(name+'-result.json')
        invoke([str(solver),'solve-stages',*map(str,paths),str(target)],output/(name+'.log'))
        report=output/(name+'-comparison.json');verify(folder/'stage-results.json',target,report)
        rows.append({'case':name,'exact_states_and_frames':True,'inputs_sha256':{str(p):sha(p) for p in paths}})
    with (output/'verification.json').open('x') as f:json.dump({'schema':'soil-fem-water-legacy/1','passed':True,'solver_sha256':sha(solver),'cases':rows},f,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();compare(a.solver.resolve(),a.output_dir.resolve())
