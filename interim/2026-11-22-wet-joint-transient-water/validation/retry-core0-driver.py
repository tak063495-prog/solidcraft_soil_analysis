import hashlib,json,subprocess,sys,time
from pathlib import Path
h=Path('interim/2026-11-22-wet-joint-transient-water/validation')
frozen=json.loads((h/'source-before.json').read_text())
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest for p,digest in frozen.items())
cmd=[sys.executable,'tools/run_ci_python_tests.py','--group','python-core-0','--output',str(h/'python-core-0-retry')]
start=time.monotonic()
with (h/'python-core-0-retry.console.log').open('x') as f:
 child=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT);code=child.wait()
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest for p,digest in frozen.items())
(h/'python-core-0-retry-process.json').write_text(json.dumps(dict(name='python-core-0',command=cmd,pid=child.pid,exit_code=code,elapsed_seconds=time.monotonic()-start,reason='Original archived evidence was not yet materialized; downloaded and SHA-verified the unchanged Git blobs. Entire 274-test group rerun without test/code/tolerance changes.',original_failed_output='python-core-0',source_hashes_unchanged=True),indent=2)+'\n')
assert code==0
