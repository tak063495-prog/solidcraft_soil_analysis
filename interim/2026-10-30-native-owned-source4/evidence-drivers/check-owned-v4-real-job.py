"""Compile and run a fresh real public-owned workflow after serialized big jobs."""
from pathlib import Path
import datetime, hashlib, json, os, shutil, subprocess, sys, time

p = Path(__file__).resolve().parent
root = p.parents[1]
out = p / 'native-owned-v4-real-job-v1'
out.mkdir()
read = lambda q: json.loads(q.read_text('utf-8-sig'))

def sha(q):
    with q.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def dump(q, value):
    with q.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False)

def pause_check():
    if (p / 'pause-requested').exists():
        raise SystemExit(2)

def gate(q):
    while not q.exists():
        pause_check()
        time.sleep(5)
    assert read(q)['passed'], str(q)

pause_check()
dump(out / 'driver-child.json', dict(pid=os.getpid(), executable=sys._base_executable,
    script=str(Path(__file__).resolve()), created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump(out / 'protocol.json', dict(prepared=True, actual_real_FEM_worker_not_stub=True,
    intended_mesh_target_m=2., element_order=2, quadrature_order=3, checkpoint_interval=4,
    full_prior_large_audits_and_surrounding16088_and_native_CI_serialized_before_new_solve=True,
    public_start_status_complete_open_configure_query_required=True,
    exact_job_safe_pause_registration_required=True, no_deadline_kill_of_numerical_child=True,
    large14909_GUI_and_general_convergence_not_established_by_this_small_case=True, p6_complete=False))
gate(p / 'native-owned-v4-public-fixed-v1/terminal.json')
gate(p / 'native-owned-v4-quality-v1/terminal.json')
host_parent = p / 'native-owned-v4-real-host-v1'
assert read(host_parent / 'terminal.json')['passed']
host = host_parent / 'runtime'
native = p / 'native-owned-v4'
library_path = p / 'native-owned-v4-public-fixed-v1/library-binding.json'
library = read(library_path)
assert library['passed'] and library['selected_by_actual_engine_dependency']
engine, serde = Path(library['engine']), Path(library['serde'])
assert sha(engine) == library['engine_sha256'] and sha(serde) == library['serde_sha256']
native_sources = read(p / 'native-owned-v4-quality-v1/source-manifest.json')['files_sha256']
host_sources = read(host_parent / 'source-manifest.json')
document = p / 'actual14909-blocked-v2/workflow/conformal/source.solidcraft.json'
protocol = read(host_parent / 'protocol.json')
assert sha(document) == protocol['source_document_sha256']
solver = host / 'target/release/soil-fem-3d.exe'
assert sha(solver) == protocol['solver_sha256']
bindings = {str(native / name): h for name, h in native_sources.items()}
bindings.update({str(host / name): h for name, h in host_sources.items()})
bindings.update({str(q): sha(q) for q in [document, solver, library_path, engine, serde,
    p / 'owned-v4-real-job-probe.rs', Path(__file__)]})
dump(out / 'source-binding.json', bindings)
env = dict(os.environ, PYTHONUTF8='1', PYTHONDONTWRITEBYTECODE='1',
    OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1',
    SOLIDCRAFT_SOIL_FEM_CHECKPOINT_INTERVAL='4')
for name in ('PYTHONPATH', 'PYTHONHOME'):
    env.pop(name, None)
records = []

def call(label, command):
    pause_check()
    assert all(sha(Path(q)) == h for q, h in bindings.items())
    start = time.monotonic()
    with (out / (label+'.log')).open('x', encoding='utf-8') as log:
        child = subprocess.Popen(list(map(str, command)), cwd=native, env=env, stdout=log,
            stderr=subprocess.STDOUT, creationflags=subprocess.CREATE_NO_WINDOW | subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out / (label+'-child.json'), dict(pid=child.pid, parent_pid=os.getpid(), executable=str(command[0]),
            arguments=list(map(str, command[1:])), created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code = child.wait()
    assert all(sha(Path(q)) == h for q, h in bindings.items())
    row = dict(step=label, exit_code=code, elapsed_seconds=time.monotonic()-start, bound_sources_unchanged=True)
    records.append(row)
    dump(out / (label+'-terminal.json'), row)
    print(json.dumps(row), flush=True)
    if code == 3:
        assert read(out / 'run/stopped.json')['cooperative_stop']
        raise SystemExit(3)
    assert code == 0, label

binary = out / 'owned-v4-real-job-probe.exe'
call('compile', ['C:/Users/link_/.cargo/bin/rustc.exe', '--edition', '2024',
    '-L', 'dependency='+str(serde.parent), '--extern', 'solidcraft_engine='+str(engine),
    '--extern', 'serde_json='+str(serde), p / 'owned-v4-real-job-probe.rs', '-o', binary])
for item in ['audit-memory-real14909-v2/terminal.json', 'audit-memory-real14909-v2/worker-passed.json',
             'geometry-audit-real8061-v1/terminal.json',
             'global-volume-sensitivity-study-v1/terminal.json', 'native-owned-v4-main-ci-v1/terminal.json']:
    gate(p / item)
policy_path = root / 'runs/goal-progress-20261027-v1/active-pause-policy.json'
policy = read(policy_path)
guard = Path(policy['script'])
assert sha(guard) == policy['script_sha256']
job = out / 'run/owned-case/job/workflow/analysis/solver-job'
relative = job.relative_to(root).as_posix()
assert "jobs.append(root/'"+relative+"')" in guard.read_text('utf-8')
while True:
    pause_check()
    memory = json.loads(subprocess.check_output(['pwsh', '-NoProfile', '-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
        encoding='utf-8-sig', creationflags=subprocess.CREATE_NO_WINDOW))
    if min(memory['FreePhysicalMemory'], memory['FreeVirtualMemory'])*1024 >= 20*(1 << 30) and shutil.disk_usage(host).free >= 32*(1 << 30):
        break
    time.sleep(5)
dump(out / 'resource-preflight.json', dict(passed=True, minimum_each_free_bytes=20*(1 << 30),
    free_physical_bytes=memory['FreePhysicalMemory']*1024, free_virtual_bytes=memory['FreeVirtualMemory']*1024,
    free_disk_bytes=shutil.disk_usage(host).free, serialized_prior_big_phases_completed=True,
    exact_solver_job_registered_before_launch=True, solver_job=str(job), guard_sha256=sha(guard)))
call('real-owned-workflow', [binary, document, host, out / 'run', p / 'pause-requested'])
receipt = read(out / 'run/receipt.json')
assert receipt['passed'] and receipt['view4_all5_frames']
manifest = read(out / 'run/owned-case/job/workflow/view/view-manifest.json')
assert manifest['schema'] == 'soil-fem-solidcraft-view/4' and manifest['complete']
assert manifest['semantics']['raw_point_storage'] == 'authoritative_source'
assert manifest['semantics']['source_point_descriptor_schema'] == 'soil-fem-view-source-points/1'
audit = read(out / 'run/owned-case/job/workflow/view-verification.json')
assert audit['passed']
dump(out / 'terminal.json', dict(passed=True, records=records, real_public_owned_MCC_and_audits=True,
    all5_source4_public_frames_open_configure_query=True, original_bound_sources_unchanged=True,
    small_workflow_only=True, actual_GUI_and_large_case_and_spatial_convergence_separate=True, p6_complete=False))
print('PASS real public-owned small MCC workflow; all P6 remains separate', flush=True)
