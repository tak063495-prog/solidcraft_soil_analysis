"""Run every existing evidence shard, two independent processes at a time."""
import datetime,json,subprocess,sys,time
from pathlib import Path
root=Path('/workspace/work/a01-source');out=Path('/workspace/work/a01-full-validation')
queue=['python-mean-evidence','python-core-0','python-core-1','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause']
active={};results=[];started=time.monotonic()
def save():
 (out/'python-driver-progress.json').write_text(json.dumps(dict(completed=results,running=list(active),pending=queue,elapsed_seconds=time.monotonic()-started),indent=2)+'\n')
while queue or active:
 while queue and len(active)<2:
  group=queue.pop(0);log=(out/(group+'-driver.log')).open('x',encoding='utf-8')
  command=[sys.executable,str(root/'tools/run_ci_python_evidence_shards.py'),'--group',group,'--output',str(out/group)]
  child=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT)
  active[group]=(child,log,time.monotonic(),command);print('Started',group,child.pid,flush=True);save()
 for group,(child,log,t0,command) in list(active.items()):
  code=child.poll()
  if code is None:continue
  log.close();results.append(dict(group=group,pid=child.pid,exit_code=code,elapsed_seconds=time.monotonic()-t0,command=command))
  del active[group];print('Finished',group,'exit',code,flush=True);save()
 time.sleep(1)
save();(out/'python-driver-results.json').write_text(json.dumps(results,indent=2)+'\n')
sys.exit(int(any(r['exit_code'] for r in results)))
