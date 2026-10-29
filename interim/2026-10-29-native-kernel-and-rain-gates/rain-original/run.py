"""Frozen hydraulic comparison, with independent reference checks first."""
from pathlib import Path
import copy
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile

p = Path(__file__).resolve().parent
root = p.parents[1]
read = lambda q: json.loads(q.read_text('utf-8-sig'))
sha = lambda q: hashlib.sha256(q.read_bytes()).hexdigest()
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONUTF8='1', PYTHONDONTWRITEBYTECODE='1')


def write(path, value):
    with path.open('x', encoding='utf-8') as f:
        json.dump(value, f, allow_nan=False, indent=2)


def call(label, args, folder=p):
    if (p/'pause-requested').exists():
        raise SystemExit(2)
    started = time.monotonic()
    command = list(map(str, args))
    with (folder/(label+'.log')).open('x', encoding='utf-8') as log:
        child = subprocess.Popen(command, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW | subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(folder/(label+'-child.json'), dict(pid=child.pid, parent_pid=os.getpid(), executable=command[0],
            arguments=command[1:], created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code = child.wait()
    receipt = dict(step=label, exit_code=code, elapsed_seconds=time.monotonic()-started)
    write(folder/(label+'-terminal.json'), receipt)
    print(json.dumps(receipt), flush=True)
    assert code == 0, receipt


write(p/'driver-child.json', dict(pid=os.getpid(), executable=sys._base_executable,
    script=str(Path(__file__).resolve()), created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
# Every input and procedure is fixed before any independent/native comparison.
sys.path.insert(0, str(root/'tools'))
from freeze_sources import freeze
freeze(p)
runtime = p/'runtime'
runtime.mkdir()
with zipfile.ZipFile(p/'source-snapshot.zip') as archive:
    archive.extractall(runtime)
assert len(read(p/'source-snapshot.json')) == 611
sys.path.insert(0, str(runtime/'tools'))
procedures = {q.name: sha(q) for q in p.glob('*.py')}
write(p/'procedures.json', procedures)
call('reference-tests', [sys.executable, '-X', 'utf8', '-m', 'unittest', 'discover', '-s', p, '-p', 'test_reference.py', '-v'])
assert 'Ran 3 tests' in (p/'reference-tests.log').read_text('utf-8')
solver = p/'solver/soil-fem-3d.exe'
solver.parent.mkdir()
origin = root/'target/release/soil-fem-3d.exe'
assert sha(origin) == 'd55ffe42c3b61e0891fec513489c9f5881830733ac9e003a01a0cf69f495905d'
shutil.copy2(origin, solver)
base = read(runtime/'examples/richards-rain-pause.json')
write(p/'reference-case.json', base)
write(p/'scope.json', dict(schema='soil-fem-local-rain-study/1', source_files=611,
    solver_sha256=sha(solver), rain_patch_m=[[.25,.75],[.25,.75]], rain_area_m2=.25,
    geometry='1m cube, bottom H=0, sides no-flow, top atmosphere with a central square rain patch',
    physical_periods=base['analysis']['periods'], native_subdivisions=[4,8,16], native_time_step_s=[100,50,25],
    independent_subdivisions=[8,16,32], independent_time_step_s=[100,50,25],
    predeclared_gates=dict(finest_head_rms_m=.003, finest_storage_error_m3=.0001,
        successive_space_and_time_difference_ratio=.8, independent_global_mass_error_m3=1e-9),
    mechanics=False, surface_storage_or_runoff_dynamics=False, current_p6_complete=False,
    limits='Homogeneous/isotropic transient water-only benchmark; no arbitrary geology, plastic settlement or rainfall stability certification.'))
for n, dt in [(8,100),(16,100),(32,100),(16,50),(16,25)]:
    call(f'fv-{n}-dt{dt}', [sys.executable,'-X','utf8',p/'reference.py',p/'reference-case.json',p/f'fv-{n}-dt{dt}.json','--n',n,'--dt',dt])

from richards_dam_mesh import make
from collections import Counter
for n, dt in [(4,100),(8,100),(16,100),(8,50),(8,25)]:
    if (p/'pause-requested').exists():
        raise SystemExit(2)
    folder = p/f'rt0-{n}-dt{dt}'
    folder.mkdir()
    mesh = make(n)
    faces = Counter(tuple(sorted(t[i] for i in face)) for t in mesh['tetrahedra'] for face in [(1,2,3),(0,3,2),(0,1,3),(0,2,1)])
    mesh['face_sets'] = dict(bottom=[], rain=[], dry=[])
    for face, count in faces.items():
        if count != 1:
            continue
        coords = [mesh['nodes'][i] for i in face]
        if all(q[2] == 0. for q in coords):
            mesh['face_sets']['bottom'].append(list(face))
        elif all(q[2] == 1. for q in coords):
            patch = all(.25 <= q[0] <= .75 and .25 <= q[1] <= .75 for q in coords)
            mesh['face_sets']['rain' if patch else 'dry'].append(list(face))
    case = copy.deepcopy(base)
    case['settings'] = dict(initial_time_step_s=float(dt),maximum_time_step_s=float(dt),
        absolute_flow_tolerance_m3_s=1e-17,store_integration_points=False)
    for period in case['analysis']['periods']:
        top = period['boundaries'].pop()
        assert top['face_set'] == 'top'
        top['face_set'] = 'rain'
        period['boundaries'].extend([top,dict(face_set='dry',condition=dict(kind='atmospheric_no_ponding',rainfall_vertical_m_s=0.))])
    write(folder/'mesh.json', mesh)
    write(folder/'case.json', case)
    write(folder/'inputs.json', dict(mesh_sha256=sha(folder/'mesh.json'),case_sha256=sha(folder/'case.json'),solver_sha256=sha(solver)))
    call('solve', [solver,'solve-richards',folder/'mesh.json',folder/'case.json',folder/'result.json'], folder)
    run = read(folder/'result.json')
    assert run['complete'] and run['solver_binary_sha256'] == sha(solver)
    call('discrete-audit', [sys.executable,'-X','utf8',runtime/'tools/verify_richards.py',folder/'mesh.json',folder/'result.json',folder/'audit.json'],folder)
    write(folder/'terminal.json', dict(passed=True,all_periods_complete=True,independent_discrete_audit=True,
        independent_pde_and_refinement_comparison_pending=True))
assert procedures == {name: sha(p/name) for name in procedures}
write(p/'raw-terminal.json', dict(passed=True,reference_tests=3,independent_runs=5,native_runs=5,
    comparisons_and_checkpoint_controls_pending=True, p6_complete=False))
print('PASS raw local-rain solves and discrete audits; final comparisons pending',flush=True)
