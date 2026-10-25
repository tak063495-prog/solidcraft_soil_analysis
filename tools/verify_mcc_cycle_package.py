"""Reassemble the completed dry MCC cycle and verify release and VTU data."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from high_order_io import stages
from verify_high_order import verify_stages


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(out):
    root=Path(__file__).resolve().parents[1];run=root/'runs/mcc-rate-general-development-20261024-v1'
    source=root/'runs/geostatic-test-20261024-003843-2c98f3fd/slope-pit/analysis/ground'
    solver=root/'target/release/soil-fem-3d.exe'
    inputs=[source/'mesh.json',source/'ground.json',run/'slope-cycle-global-refined-case.json',run/'slope-cycle-rust-global-refined.json',
            solver,Path(__file__),root/'tools/high_order_io.py',root/'tools/verify_high_order.py',root/'tools/verify_geostatic.py',root/'tools/fem_reference.py']
    fingerprints={str(p):sha(p) for p in inputs}
    out.mkdir(exist_ok=False)
    for p,name in zip(inputs[:4],('mesh.json','ground.json','stages.json','result.json')):shutil.copyfile(p,out/name)
    stages(out/'mesh.json',out/'result.json',out/'raw-vtu')
    verify_stages(out/'mesh.json',out/'ground.json',out/'result.json',out/'raw-vtu',solver)
    with (out/'stage-verification.json').open() as f:checked=json.load(f)
    assert len(checked['frames'])==4 and checked['raw_integration_point_roundtrip']=='passed'
    assert fingerprints=={str(p):sha(p) for p in inputs},'source/input/solver changed during audit'
    record=dict(passed=True,goal_complete=False,p6_complete=False,
        scope='Independent seed/internal/physical/holding-force assembly, dry MCC cycle and raw VTU readback; independent nonlinear solution is a separate gate',
        frames=checked['frames'],excavation_release_and_removed_history_checked=True,input_and_source_sha256=fingerprints,
        evidence_sha256={p.relative_to(out).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'verification.json').open('x') as f:json.dump(record,f,indent=2)
    print('Dry MCC cycle force, excavation release, monitors and raw VTU audit passed',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();verify(a.output.resolve())
