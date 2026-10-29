from pathlib import Path
import datetime, hashlib, json, os, subprocess, sys, time
p=Path(__file__).resolve().parent
read=lambda q:json.loads(q.read_text('utf-8-sig'))
out=p/'public14909-reader-v3'
assert read(p/'actual14909-native-package-repair-v1/terminal.json')['passed']
prep=read(out/'preparation.json')
assert prep['prepared']
script=p/'read-real14909-v3.py'
with script.open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==prep['script_sha256']
if (p/'pause-requested').exists():raise SystemExit(2)
with (out/'driver-child.json').open('x',encoding='utf-8') as f:
    json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f,indent=2)
env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
start=time.monotonic()
with (out/'data-phase.log').open('x',encoding='utf-8') as log:
    child=subprocess.Popen([sys.executable,'-I','-X','utf8',str(script)],cwd=p.parents[1],env=env,
        stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    with (out/'data-phase-child.json').open('x',encoding='utf-8') as f:
        json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=sys.executable,
            arguments=['-I','-X','utf8',str(script)],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f,indent=2)
    code=child.wait()
with (out/'process-terminal.json').open('x',encoding='utf-8') as f:
    json.dump(dict(passed=code==0,exit_code=code,elapsed_seconds=time.monotonic()-start,p6_complete=False),f,indent=2)
assert code==0
assert read(out/'terminal.json')['passed']
print('PASS actual14909 original public Native reader after verified canonical audit placement',flush=True)
