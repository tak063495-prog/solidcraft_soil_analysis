"""Phase1 only: preserve an accepted seed, then recompute remaining water time."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];folder=p/'rt0-24-dt100'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    h=hashlib.sha256()
    with q.open('rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()
def write(q,v):
    with q.open('x') as f:json.dump(v,f,indent=2)
write(p/'native-driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
protocol=read(p/'protocol.json');solver=p/'solver/soil-fem-3d.exe';assert sha(solver)==protocol['solver_sha256']
binding={n:sha(folder/n) for n in ['mesh.json','case.json']}
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
def call(label,args):
    if (p/'pause-requested').exists():raise SystemExit(2)
    start=time.monotonic();command=list(map(str,args))
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(folder/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
    r=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-start);write(folder/(label+'-terminal.json'),r);print(json.dumps(r),flush=True);assert code==0
while True:
    if (p/'pause-requested').exists():raise SystemExit(2)
    free=int(subprocess.check_output(['pwsh','-NoProfile','-Command','(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory'],encoding='utf-8-sig').strip())*1024
    if free>=24*(1<<30):break
    time.sleep(5)
write(folder/'resource-preflight.json',dict(free_physical_bytes=free,required_free_bytes=24*(1<<30),no_solver_limits_or_tolerances_raised=True))
call('seed-checkpoint',[solver,'checkpoint-richards',folder/'mesh.json',folder/'case.json','1',folder/'seed-checkpoint.json'])
seed_sha=sha(folder/'seed-checkpoint.json')
write(folder/'seed-identity.json',dict(checkpoint_sha256=seed_sha,solver_sha256=sha(solver),input_sha256=binding,accepted_steps=1))
call('resume-transient',[solver,'resume-richards',folder/'mesh.json',folder/'case.json',folder/'seed-checkpoint.json',folder/'result.json'])
assert sha(folder/'seed-checkpoint.json')==seed_sha and binding=={n:sha(folder/n) for n in binding}
write(folder/'raw-result-terminal.json',dict(passed=True,solver_sha256=sha(solver),result_sha256=sha(folder/'result.json'),result_bytes=(folder/'result.json').stat().st_size,
    first_accepted_seed_preserved=True,independent_discrete_audit_and_refinement_acceptance_pending=True,p6_complete=False))
print('Phase1 native24 raw complete; independent audit remains required',flush=True)
