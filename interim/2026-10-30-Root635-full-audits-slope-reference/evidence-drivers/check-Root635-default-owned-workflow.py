"""Exercise the adopted, unmodified Root through the actual owned Native worker."""
from pathlib import Path
import datetime, hashlib, json, os, shutil, subprocess, sys, time

p = Path(__file__).resolve().parent
root = p.parents[1]
out = p / 'Root635-default-owned-workflow-v1'
read = lambda q: json.loads(q.read_text('utf-8-sig'))

def sha(q):
    with q.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def dump(q, value):
    with q.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False)

def active():
    if (p / 'pause-requested').exists():
        raise SystemExit(2)

def gate(q):
    while not q.exists():
        active()
        time.sleep(5)
    assert read(q)['passed'], str(q)

active()
prep = read(out / 'preparation.json')
assert prep['prepared'] and prep['actual_Root_files635_unchanged']
binding = read(out / 'source-binding.json')
sources = binding['Root_files_sha256']
binary = Path(prep['binary'])
document = Path(prep['document'])
assert Path(prep['host']) == root and len(sources) == 635
assert sha(binary) == prep['binary_sha256']
assert sha(document) == binding['original_document_sha256']
dump(out / 'workflow-driver-child.json', dict(pid=os.getpid(), executable=sys._base_executable,
    script=str(Path(__file__).resolve()), created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump(out / 'workflow-protocol.json', dict(actual_adopted_Root=True, source_files=635,
    no_source_point_forcing_or_private_test_worker=True, default_mesh_target_m=2.0,
    element_order=2, quadrature_order=3, checkpoint_interval=4,
    full_actual_solve_and_original_worker_audits_required=True,
    five_public_frames_required=True, source_files_frozen_until_terminal=True,
    small_case_not_general_convergence_or_new_Native_GUI_acceptance=True, p6_complete=False))

# Registration is supplied only after the second increment-control driver has
# consumed its pinned old registration. Waiting creates no numerical job.
gate(out / 'pause-registration.json')
gate(p / 'native-owned-v4-main-ci-v1/terminal.json')
policy = read(root / 'runs/goal-progress-20261027-v1/active-pause-policy.json')
guard = Path(policy['script'])
assert sha(guard) == policy['script_sha256']
job = out / 'run/owned-case/job/workflow/analysis/solver-job'
relative = job.relative_to(root).as_posix()
registration = read(out / 'pause-registration.json')
assert registration['job'] == relative
assert "jobs.append(root/'" + relative + "')" in guard.read_text('utf-8')

while True:
    active()
    memory = json.loads(subprocess.check_output(['pwsh', '-NoProfile', '-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
        encoding='utf-8-sig', creationflags=subprocess.CREATE_NO_WINDOW))
    if min(memory['FreePhysicalMemory'], memory['FreeVirtualMemory']) * 1024 >= 20 * (1 << 30) and shutil.disk_usage(root).free >= 32 * (1 << 30):
        break
    time.sleep(5)
active()
assert all(sha(root / name) == h for name, h in sources.items())
assert sha(root / 'target/release/soil-fem-3d.exe') == binding['Root_solver_sha256']
assert sha(binary) == prep['binary_sha256']
assert sha(document) == binding['original_document_sha256']
dump(out / 'resource-preflight.json', dict(passed=True,
    free_physical_bytes=memory['FreePhysicalMemory'] * 1024,
    free_virtual_bytes=memory['FreeVirtualMemory'] * 1024,
    minimum_each_free_bytes=20 * (1 << 30), free_disk_bytes=shutil.disk_usage(root).free,
    exact_job_registered_before_launch=True, job=str(job), guard_sha256=sha(guard)))
env = dict(os.environ, PYTHONUTF8='1', PYTHONDONTWRITEBYTECODE='1',
    OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1', RAYON_NUM_THREADS='1',
    SOLIDCRAFT_SOIL_FEM_CHECKPOINT_INTERVAL='4')
for key in ('PYTHONPATH', 'PYTHONHOME'):
    env.pop(key, None)
command = [str(binary), str(document), str(root), str(out / 'run'), str(p / 'pause-requested')]
start = time.monotonic()
with (out / 'real-owned-workflow.log').open('x', encoding='utf-8') as log:
    child = subprocess.Popen(command, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT,
        creationflags=subprocess.CREATE_NO_WINDOW | subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    dump(out / 'real-owned-workflow-child.json', dict(pid=child.pid, parent_pid=os.getpid(),
        executable=command[0], arguments=command[1:],
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    code = child.wait()
unchanged = all(sha(root / name) == h for name, h in sources.items())
assert unchanged
dump(out / 'real-owned-workflow-terminal.json', dict(exit_code=code,
    elapsed_seconds=time.monotonic()-start, all635_bound_Root_files_unchanged=unchanged))
if code == 3:
    assert read(out / 'run/stopped.json')['cooperative_stop']
    raise SystemExit(3)
assert code == 0
receipt = read(out / 'run/receipt.json')
assert receipt['passed'] and receipt['default_view_all5_frames']
assert receipt['actual_adopted_Root_unmodified'] and len(receipt['frames']) == 5
assert receipt['original_projected_CAD_input_sha256'] == prep['input_sha256']
status = read(out / 'run/terminal-status.json')
assert status['state'] == 'succeeded' and status['result_usable']
manifest = read(Path(status['result']))
assert manifest['schema'] in ('soil-fem-solidcraft-view/3', 'soil-fem-solidcraft-view/4')
assert manifest['complete'] and len(manifest['frames']) == 5
audit = read(out / 'run/owned-case/job/workflow/view-verification.json')
assert audit['passed']
dump(out / 'terminal.json', dict(passed=True, actual_Root_files=635,
    all_bound_Root_files_unchanged=True, actual_public_owned_MCC_solve_and_original_audits=True,
    actual_default_view_schema=manifest['schema'], default_format_not_forced=True,
    all5_public_frames_open_configure_query=True, actual_Root_worker_and_numerics_unmodified=True,
    actual_new_GUI_and_spatial_convergence_and_new_PC_and_allP6_separate=True, p6_complete=False))
print('PASS actual adopted Root default public owned MCC workflow', flush=True)
