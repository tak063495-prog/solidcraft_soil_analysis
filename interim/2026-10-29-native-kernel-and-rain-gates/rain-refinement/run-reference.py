"""Independent FV refinement; no native head/stiffness/result is read."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1]
def write(q,v):
    with q.open('x') as f:json.dump(v,f,indent=2)
write(p/'reference-driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
records=[]
for n,dt in [(64,100),(32,50),(32,25)]:
    if (p/'pause-requested').exists():raise SystemExit(2)
    label=f'fv-{n}-dt{dt}';command=list(map(str,[sys.executable,'-X','utf8',p/'reference.py',p/'reference-case.json',p/(label+'.json'),'--n',n,'--dt',dt]));start=time.monotonic()
    with (p/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(p/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-start);write(p/(label+'-terminal.json'),row);records.append(row);print(json.dumps(row),flush=True);assert code==0
write(p/'reference-terminal.json',dict(passed=True,independent_pde_runs=3,records=records,comparison_and_mechanics_not_accepted=True,p6_complete=False))
print('Phase1 independent FV refinement complete; cross-solver/refinement gates remain required',flush=True)
