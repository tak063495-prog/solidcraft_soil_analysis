"""Retain a rejected elastic formulation and replay its inputs with the fix."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from verify_checkpoint import invoke


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(solver,old,output):
    output.mkdir(parents=True,exist_ok=False)
    evidence=read(old/'counterexample.json')
    assert evidence['counterexample_confirmed'] and not evidence['formulation_accepted']
    for path,digest in evidence['sha256'].items():assert sha(old/path)==digest
    assert sha(old/'solver/solver.exe')==evidence['solver_sha256']
    shutil.copytree(old,output/'rejected-formulation')
    rows=[]
    for label in ('1e-04','1e-05'):
        inp=output/'rejected-formulation'/(label+'-input.json');out=output/(label+'-corrected.json')
        invoke([str(solver),'hardening-path',str(inp),str(out)],output/(label+'.log'))
        d=read(out);seed=read(inp)['initial_state'];last=d['states'][-1]
        error=float(np.max(np.abs(np.asarray(last['effective_stress_compression_positive_kpa'])-seed['effective_stress_compression_positive_kpa'])))
        assert d['plastic_substeps']==0 and error<1e-9
        assert abs(last['void_ratio']-seed['void_ratio'])<1e-12
        assert last['cap_pressure_kpa']==seed['cap_pressure_kpa']
        assert last['plastic_strain_compression_positive']==seed['plastic_strain_compression_positive']
        rows.append({'maximum_substep':read(inp)['maximum_substep'],'maximum_residual_stress_kpa':error,'plastic_substeps':0})
    report={'schema':'soil-fem-hardening-elastic-cycle-verification/1','passed':True,'solver_sha256':sha(solver),
        'rejected_solver_sha256':evidence['solver_sha256'],'corrected_cycles':rows,
        'rejected_formulation_evidence_sha256':sha(output/'rejected-formulation/counterexample.json')}
    with (output/'verification.json').open('x') as f:json.dump(report,f,indent=2)
    print('Corrected nonproportional elastic cycles return to initial state:',rows)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--rejected',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.rejected.resolve(),a.output_dir.resolve())
