"""Bounded independent HTTP/1.1 fetches; never fabricate a missing artifact."""
from pathlib import Path
import subprocess,os,json,time,concurrent.futures
root=Path.cwd();env=dict(os.environ,GIT_NO_LAZY_FETCH='1');base='11df919ad053fc742fa40cb986953d8f53de28db';rows=[];oid_paths={}
for line in subprocess.check_output(['git','ls-tree','-r',base],text=True,env=env).splitlines():
 meta,name=line.split('\t');mode,kind,oid=meta.split()
 if kind!='blob' or name.startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part'):continue
 oid_paths.setdefault(oid,[]).append(name)
oids=sorted(oid_paths)
def check():
 p=subprocess.run(['git','cat-file','--batch-check=%(objectname) %(objecttype)'],input='\n'.join(oids)+'\n',text=True,capture_output=True,env=env,check=True)
 return [line.split()[0] for line in p.stdout.splitlines() if line.endswith(' missing')]
missing=check();print('remaining required Git blobs',len(missing),flush=True)
logs=root/'work/a02follow-small-fetches-v2';logs.mkdir(exist_ok=True)
def fetch(oid):
 attempts=[]
 for index in range(3):
  cmd=['git','-c','http.version=HTTP/1.1','-c','credential.helper=!gh auth git-credential','-c','fetch.negotiationAlgorithm=noop','-c','gc.auto=0','-c','maintenance.auto=false','fetch','--no-tags','--no-write-fetch-head','--filter=blob:none','origin',oid]
  start=time.monotonic()
  with (logs/(oid+f'-{index}.log')).open('x') as log:
   p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT);timed=False
   try:code=p.wait(timeout=300)
   except subprocess.TimeoutExpired:timed=True;p.kill();code=p.wait()
  attempts.append(dict(command=cmd,pid=p.pid,exit_code=code,timed_out=timed,elapsed_seconds=time.monotonic()-start))
  if code==0:break
 result=dict(oid=oid,paths=oid_paths[oid],attempts=attempts,success=code==0);(logs/(oid+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(oid,code,flush=True);return result
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(fetch,missing))
remaining=check();v=dict(initial_missing=len(missing),attempts=rows,remaining=remaining,all_required_objects_available=not remaining,excluded_only_existing_postfold=True,working_files_not_replaced_with_synthetic_values=True)
(root/'work/a02follow-small-fetch-report.json').write_text(json.dumps(v,indent=2)+'\n');assert not remaining
