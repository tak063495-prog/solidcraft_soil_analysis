"""Compare old/new local-excavation solvers on the same freshly generated CAD mesh."""
import argparse
import hashlib
import json
from pathlib import Path
from verify_checkpoint import invoke
from verify_mcc_regression import verify


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def compare(baseline,current,output,comparison):
    output.mkdir(parents=True,exist_ok=False)
    accepted_path=baseline/'final-validation.json'
    accepted=json.loads(accepted_path.read_text(encoding='utf-8-sig'))
    solver=baseline/'solver/soil-fem-3d.exe'
    assert accepted['passed'] and sha(solver)==accepted['solver_sha256']
    paths=[current/'ground/mesh.json',current/'ground/ground.json',current/'stages.json']
    result=output/'old-solver-same-input.json'
    invoke([str(solver),'solve-stages',*map(str,paths),str(result)],output/'solve.log')
    verify(result,current/'stage-results.json',comparison)
    proof={'schema':'soil-fem-local-same-input-regression/1','passed':True,
        'baseline_solver_sha256':sha(solver),'accepted_baseline_manifest_sha256':sha(accepted_path),
        'inputs_sha256':{str(p):sha(p) for p in paths},'exact_states_and_frames':True,
        'scope':'Same fresh CAD mesh and loads are used by both solvers. Independent CAD regeneration can change curved coordinates by rounding; no tolerance is substituted for solver-state equality.'}
    with (output/'verification.json').open('x',encoding='utf-8') as stream:json.dump(proof,stream,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,required=True);p.add_argument('--current',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--comparison',type=Path,required=True)
    a=p.parse_args();compare(a.baseline.resolve(),a.current.resolve(),a.output_dir.resolve(),a.comparison.resolve())
