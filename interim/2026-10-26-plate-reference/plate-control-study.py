"""Execute the fixed public plate/beam/soil example and all scoped audits."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

folder=Path(__file__).resolve().parent
root=folder.parents[1]
sys.path.insert(0,str(root/'tools'))
import freeze_sources

parser=argparse.ArgumentParser()
parser.add_argument('--solver',type=Path,required=True)
args=parser.parse_args()
run=folder/'control'
run.mkdir()
freeze_sources.freeze(run)
sources=json.loads((run/'source-snapshot.json').read_text())
for name in ('solver','ground'):(run/name).mkdir()
solver=run/'solver/soil-fem-3d.exe';shutil.copyfile(args.solver.resolve(),solver)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
digest=sha(solver)
(run/'solver-binding.json').write_text(json.dumps(dict(solver_sha256=digest)),encoding='utf-8')
for source,destination in [('examples/hardening-column-tet10-mesh.json','ground/mesh.json'),
    ('examples/ground-plate-beam-column.json','ground/ground.json'),
    ('examples/stages-plate-beam-column.json','stages.json')]:
    shutil.copyfile(root/source,run/destination)

def invoke(name,command):
    assert freeze_sources.manifest()==sources
    code=subprocess.call([sys.executable,'-X','utf8',str(root/'tools/run_logged.py'),str(run/(name+'.log')),*map(str,command)],cwd=root)
    if code:raise SystemExit(code)
    assert sha(solver)==digest and freeze_sources.manifest()==sources

invoke('solve',[solver,'solve-stages',run/'ground/mesh.json',run/'ground/ground.json',run/'stages.json',run/'stage-results.json'])
invoke('export',[sys.executable,'-X','utf8',root/'tools/export_stages.py',run/'ground/mesh.json',run/'stage-results.json',run/'exports'])
invoke('force-vtu-audit',[sys.executable,'-X','utf8',root/'tools/verify_stage_package.py',run/'ground/mesh.json',run/'ground/ground.json',run/'stage-results.json',run/'exports','--solver',solver])
invoke('independent',[sys.executable,'-X','utf8',root/'tools/verify_beam_coupling.py',run,run/'independent.json'])
invoke('cutback-history',[sys.executable,'-X','utf8',root/'tools/verify_beam_history.py','--solver',solver,'--run',run,'--output-dir',run/'history'])
invoke('checkpoints',[sys.executable,'-X','utf8',root/'tools/verify_construction_checkpoint.py','--solver',solver,'--run',run,'--output-dir',run/'checkpoints'])
result=json.loads((run/'stage-results.json').read_text())
assert result['complete'] and result['schema']=='soil-fem-stage-result/14'
assert len(result['frames'])==5
export=json.loads((run/'exports/export-verification.json').read_text())
assert export['plate_sample_roundtrip']=='passed'
report=dict(schema='soil-fem-plate-beam-control-verification/1',passed=True,
    scope='Explicit elastic MC branch, two flat plates and a shared-node beam, staged installation/self-weight/loading/unloading/removal',
    solver_sha256=digest,source_files=len(sources),all_3d_fields_independently_checked=True,
    raw_plate_and_soil_vtu_roundtrip=True,full_non_time_history_exact=True,
    spatial_convergence_accepted=False,actual_cad_wall_accepted=False,goal_complete=False,p6_complete=False,
    evidence_sha256={name:sha(run/name) for name in ['source-snapshot.json','stage-results.json','independent.json',
        'stage-verification.json','exports/export-verification.json','history/verification.json','checkpoints/verification.json']})
(run/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report),flush=True)
