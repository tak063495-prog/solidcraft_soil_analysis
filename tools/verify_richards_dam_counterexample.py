"""Replay the exact finite-width failure with sealed old and current engines."""
import argparse
import json
from pathlib import Path
import shutil
from verify_acceptance_packet import verify as packet,sha
from verify_checkpoint import invoke
from richards_dam_audit import audit


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def verify(solver,output):
    root=Path(__file__).resolve().parents[1];old=root/'runs/richards-finite-width-counterexample-20261024';parent=root/'runs/richards-acceptance-20261024-104752-1c503b7b'
    proof=packet(parent);assert sha(old/'before.exe')==proof['solver_sha256']
    evidence=read(old/'evidence.json')
    for n,d in evidence['sha256'].items():assert sha(old/n)==d,n
    output.mkdir(parents=True,exist_ok=False)
    for n in ('mesh.json','case.json','before.exe','evidence.json'):shutil.copyfile(old/n,output/n)
    invoke([str(output/'before.exe'),'solve-richards',str(output/'mesh.json'),str(output/'case.json'),str(output/'before.json')],output/'before.log',expected=1)
    invoke([str(solver),'solve-richards',str(output/'mesh.json'),str(output/'case.json'),str(output/'after.json')],output/'after.log')
    before=read(output/'before.json');after=read(output/'after.json')
    assert not before['complete'] and 'line search failed' in before['failure_reason'] and len(before['frames'])==1
    assert before['last_committed']==before['frames'][0]['state']==read(old/'before.json')['last_committed']
    assert before['mesh_sha256']==after['mesh_sha256'] and before['case_sha256']==after['case_sha256']
    assert after['complete'] and after['solver_binary_sha256']==sha(solver)
    checked=audit(read(output/'mesh.json'),after)
    report=dict(schema='soil-fem-richards-dam-counterexample-verification/1',passed=True,solver_sha256=sha(solver),old_parent=proof,old_parent_manifest=str((parent/'final-validation.json').relative_to(root)),exact_mesh_case_hashes_match=True,failed_old_state_not_committed=True,independent_fixed_local_audit=checked)
    with (output/'verification.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
