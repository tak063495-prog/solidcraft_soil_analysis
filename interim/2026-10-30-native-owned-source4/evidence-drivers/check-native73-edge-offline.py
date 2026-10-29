"""Supplement the skipped browser gate using installed Edge and a fresh profile."""
from pathlib import Path
import datetime, hashlib, json, os, subprocess, sys, time

p = Path(__file__).resolve().parent
root = p.parents[1]
native = root / 'runs/native-typed-source-integration-20261029-v1/native'
out = p / 'native73-edge-offline-v1'
out.mkdir()
read = lambda q: json.loads(q.read_text('utf-8-sig'))

def sha(q):
    with q.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def dump(q, value):
    with q.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False)

assert read(p / 'capacity-desktop-and-main-ci-v1/terminal.json')['passed']
assert not (p / 'pause-requested').exists()
browser = Path(r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')
assert browser.is_file()
sources = read(p / 'native-catalog-capacity-repair-v1/source-manifest.json')['files_sha256']
binding = {str(native / name): h for name, h in sources.items()}
binding[str(browser)] = sha(browser)
reference = read(p / 'public14909-reader-v2/compiled-probe-binding.json')
binding[reference['engine']] = reference['engine_sha256']
binding[reference['serde']] = reference['serde_sha256']
assert all(sha(Path(q)) == h for q, h in binding.items())
dump(out / 'source-binding.json', binding)
dump(out / 'driver-child.json', dict(pid=os.getpid(), executable=sys._base_executable,
    script=str(Path(__file__).resolve()), created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
temporary = out / 'tmp'
temporary.mkdir()
assert temporary.resolve().is_relative_to(p.resolve())
env = dict(os.environ, CARGO_TARGET_DIR=str(native/'target'), CARGO_BUILD_JOBS='1', CHROME=str(browser),
    TEMP=str(temporary), TMP=str(temporary), GIT_OPTIONAL_LOCKS='0')
command = ['C:/Users/link_/.cargo/bin/cargo.exe','xtask','web-offline-test']
start = time.monotonic()
with (out / 'test.log').open('x', encoding='utf-8') as log:
    child = subprocess.Popen(command, cwd=native, env=env, stdout=log, stderr=subprocess.STDOUT,
        creationflags=subprocess.CREATE_NO_WINDOW | subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    dump(out / 'test-child.json', dict(pid=child.pid, parent_pid=os.getpid(), executable=command[0],
        arguments=command[1:], created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    code = child.wait()
unchanged = all(sha(Path(q)) == h for q, h in binding.items())
log = (out / 'test.log').read_text('utf-8')
passed = code == 0 and unchanged and 'web-offline-test: passed' in log and 'skipped (no Chrome' not in log
dump(out / 'terminal.json', dict(passed=passed, exit_code=code, elapsed_seconds=time.monotonic()-start,
    browser=str(browser), browser_is_Edge_not_Google_Chrome=True,
    fresh_profile_under_owned_temporary_root=str(temporary), original_browser_profile_not_used=True,
    original_project_test_code_and_sources_and_actual_public_reader_libraries_unchanged=unchanged,
    default_stand_in_app_service_worker_offline_update_test_only=True,
    real_web_build_3D_FEM_or_new_owned_completion_GUI_not_tested_here=True,
    original_main_CI_Chrome_skip_and_nightly_counterexample_not_overwritten=True, p6_complete=False))
print(json.dumps(dict(passed=passed, exit_code=code, Edge_only_supplement=True)), flush=True)
raise SystemExit(0 if passed else 1)
