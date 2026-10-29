"""Explicit protocol-only worker: no solve, numerical receipt or real result."""
from pathlib import Path
import hashlib, json, os, sys

args = sys.argv[1:]
assert args[0] == 'run'
folder = Path(args[2])
folder.mkdir()
host = Path(__file__).resolve().parents[1]
mode = json.loads((host / 'mode.json').read_text('utf-8'))['mode']
identity = args[args.index('--expected-input-sha256')+1]
token = args[args.index('--host-token')+1]
assert '--cooperative' not in args
view = folder / 'workflow/view/view-manifest.json'
view.parent.mkdir(parents=True)
manifest = dict(schema='soil-fem-solidcraft-view/3' if mode == 'legacy' else
    'soil-fem-solidcraft-view/5' if mode == 'future_view' else 'soil-fem-solidcraft-view/4',
    complete=True, input_sha256=identity, solver_binary_sha256='a'*64)
if mode != 'missing_semantics':
    manifest['semantics'] = dict(source_result_schema='soil-fem-projected-mcc-stage-result/1' if mode == 'wrong_family'
        else 'soil-fem-projected-mcc-stage-result/2', raw_point_storage='averaged_points' if mode == 'wrong_storage' else 'authoritative_source',
        source_point_descriptor_schema='wrong/1' if mode == 'wrong_descriptor' else 'soil-fem-view-source-points/1')

def dump(q, value):
    with q.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2)

def sha(q):
    with q.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

dump(view, manifest)
audit = folder / 'workflow/view-verification.json'
dump(audit, dict(passed=True, input_sha256=identity, view_manifest_sha256=sha(view)))
jobid = 'explicit-protocol-test-'+mode
dump(folder / 'job.json', dict(schema='soil-fem-native-job/1', job_id=jobid, input_sha256=identity,
    host_token=token, solver_sha256='a'*64, protocol_test_only_no_numerical_solve=True))
dump(folder / 'status.json', dict(schema='soil-fem-native-job-status/1', job_id=jobid,
    input_sha256=identity, host_token=token, worker_pid=os.getpid(), state='succeeded', phase='complete',
    exit_code=0, result='workflow/view/view-manifest.json', result_sha256=sha(view), view_audit_sha256=sha(audit),
    protocol_test_only_no_numerical_solve=True))
