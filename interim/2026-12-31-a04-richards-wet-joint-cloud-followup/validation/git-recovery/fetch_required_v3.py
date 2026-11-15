from pathlib import Path
import subprocess,os,json,time,concurrent.futures
root=Path.cwd();out=Path('/workspace/work/a04follow-git-recovery');env=dict(os.environ,GIT_NO_LAZY_FETCH='1');base='d9ca6634a7c394d0188a502c11c9038a4a0b98e5';paths={}
for l in subprocess.check_output(['git','ls-tree','-r',base],text=True,env=env).splitlines():
 meta,name=l.split('\t');mode,kind,oid=meta.split()
 if kind!='blob' or name.startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part'):continue
 paths.setdefault(oid,[]).append(name)
oids=sorted(paths)
def check(ids):
 if not ids:return []
 p=subprocess.run(['git','cat-file','--batch-check=%(objectname) %(objecttype)'],input='\n'.join(ids)+'\n',text=True,capture_output=True,env=env,check=True)
 return [l.split()[0] for l in p.stdout.splitlines() if l.endswith(' missing')]
missing=check(oids);large=lambda name: any(s in name for s in ('.part','.zip','.tar','.gz','.bin','.npz'))
big=[oid for oid in missing if any(large(p) for p in paths[oid])];small=[oid for oid in missing if oid not in set(big)]
(out/'required-initial-v3.json').write_text(json.dumps(dict(base=base,total_unique_required=len(oids),initial_missing=len(missing),small=len(small),large=len(big),existing_ci_postfold_exclusion_only=True),indent=2)+'\n');print('small',len(small),'large',len(big),flush=True)
logs=out/'fetches-v3';logs.mkdir(exist_ok=True)
def fetch(job):
 ids,label=job;attempts=[]
 for retry in range(3):
  cmd=['git','-c','http.version=HTTP/1.1','-c','fetch.negotiationAlgorithm=noop','-c','gc.auto=0','-c','maintenance.auto=false','fetch','--no-tags','--no-write-fetch-head','--filter=blob:none','--stdin','origin']
  started=time.monotonic()
  with (logs/(label+'-'+str(retry)+'.log')).open('x') as f:
   p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=f,stderr=subprocess.STDOUT,text=True,env=env)
   timed=False
   try:p.communicate('\n'.join(ids)+'\n',timeout=300)
   except subprocess.TimeoutExpired:timed=True;p.kill();p.communicate()
  attempts.append(dict(command=cmd,requested=ids,pid=p.pid,exit_code=p.returncode,timed_out=timed,elapsed_seconds=time.monotonic()-started))
  if p.returncode==0:break
 result=dict(ids=ids,paths={i:paths[i] for i in ids},attempts=attempts,success=p.returncode==0)
 (logs/(label+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(label,p.returncode,flush=True);return result
rows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 rows.extend(pool.map(fetch,[(small[i:i+32],'text-'+str(i//32).zfill(3)) for i in range(0,len(small),32)]))
remaining_small=check(small);(out/'small-ready-v3.json').write_text(json.dumps(dict(all_small_available=not remaining_small,missing=remaining_small))+'\n');assert not remaining_small
with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
 rows.extend(pool.map(fetch,[([oid],'blob-'+oid) for oid in big]))
remaining=check(oids);(out/'fetch-result-v3.json').write_text(json.dumps(dict(initial_missing=len(missing),jobs=rows,remaining=remaining,all_required_objects_available=not remaining,excluded_only_existing_ci_postfold=True),indent=2)+'\n');assert not remaining
