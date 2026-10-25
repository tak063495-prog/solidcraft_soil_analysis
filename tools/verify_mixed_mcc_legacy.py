"""Execute the new development binary on accepted displacement-MCC/MC cases."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from verify_default_solver import deterministic
from verify_acceptance_packet import sha


def read(p):return json.loads(p.read_text(encoding='utf-8'))


def run(solver,out):
    root=Path(__file__).resolve().parents[1];cycle=root/'runs/mcc-tet10-cycle-acceptance-20261025-v1'
    mixed=root/'runs/mixed-mechanics-acceptance-20261025-v1';mc=root/'runs/2d-3d-column-acceptance-20261025-v1/data'
    cases=[('legacy-mcc',[cycle/'inputs/mesh.json',cycle/'inputs/ground.json',cycle/'inputs/slope-cycle-case.json'],mixed/'current-mcc/result.json'),
           ('legacy-mc',[mc/'mesh-3d.json',mc/'ground-3d.json',mc/'stages-3d.json'],mc/'result-3d.json')]
    out.mkdir(exist_ok=False);(out/'solver').mkdir();shutil.copyfile(solver,out/'solver/soil-fem-3d.exe')
    binary=out/'solver/soil-fem-3d.exe'
    paths=[solver,Path(__file__),root/'tools/verify_default_solver.py',root/'tools/run_logged.py']+[p for _,inputs,expected in cases for p in [*inputs,expected]]
    fingerprints={str(p):sha(p) for p in paths};reports=[]
    for name,inputs,expected in cases:
        output=out/(name+'.json')
        subprocess.run([sys.executable,'-X','utf8',str(root/'tools/run_logged.py'),str(out/(name+'.log')),str(binary),'solve-stages',*map(str,inputs),str(output)],check=True,cwd=root)
        a,b=read(output),read(expected)
        assert a['complete'] and a['provenance']['solver_binary_sha256']==sha(binary)
        assert json.dumps(deterministic(a),sort_keys=True,allow_nan=False)==json.dumps(deterministic(b),sort_keys=True,allow_nan=False),name
        assert all('mixed_mcc' not in f['state'] for f in a['frames'])
        reports.append(dict(name=name,all_nonruntime_fields_exact=True))
    assert fingerprints=={str(p):sha(p) for p in paths}
    report=dict(schema='soil-fem-mixed-mcc-legacy-compatibility/1',passed=True,goal_complete=False,p6_complete=False,
        scope='New development binary, accepted legacy displacement MCC/MC paths only',cases=reports,
        excluded_keys=['solver_binary_sha256','elapsed_seconds'],input_and_source_sha256=fingerprints,
        solver_sha256=sha(binary),evidence_sha256={p.relative_to(out).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'verification.json').open('x') as f:json.dump(report,f,indent=2)
    print('Legacy MCC/MC results and all non-timing diagnostics are exact.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('solver',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();run(a.solver.resolve(),a.output.resolve())
