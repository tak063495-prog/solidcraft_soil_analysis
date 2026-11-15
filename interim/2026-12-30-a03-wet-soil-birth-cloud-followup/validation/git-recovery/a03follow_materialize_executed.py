from pathlib import Path
import concurrent.futures,subprocess,hashlib,json,time
ROOT=Path.cwd();BASE='89e9dfc26cbe444199179d7faa1edb1437849309'
OUT=ROOT/'work/a03follow-materialize';OUT.mkdir(exist_ok=False)
rows=[]
for line in subprocess.check_output(['git','ls-tree','-r',BASE],text=True).splitlines():
 meta,path=line.split('\t');mode,kind,oid=meta.split()
 if kind!='blob' or path.startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part'):continue
 if not (ROOT/path).exists():rows.append((path,oid))
print('missing files',len(rows),flush=True)
def get(row):
 path,oid=row;dest=ROOT/path;dest.parent.mkdir(parents=True,exist_ok=True)
 tmp=OUT/(oid+'.download');start=time.monotonic()
 url='https://raw.githubusercontent.com/tak063495-prog/solidcraft_soil_analysis/'+BASE+'/'+path
 cmd=['curl','-fLSs','--retry','2','--max-time','300',url,'-o',str(tmp)]
 log=OUT/(oid+'.log')
 with log.open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 proof=dict(path=path,git_blob_sha1=oid,exit_code=p.returncode,seconds=time.monotonic()-start)
 if p.returncode==0:
  size=tmp.stat().st_size;h=hashlib.sha1(('blob '+str(size)+'\0').encode());s=hashlib.sha256()
  with tmp.open('rb') as f:
   for b in iter(lambda:f.read(1024*1024),b''):h.update(b);s.update(b)
  proof.update(size=size,actual_git_blob_sha1=h.hexdigest(),sha256=s.hexdigest())
  if h.hexdigest()==oid:tmp.rename(dest)
 proof['verified']=proof.get('actual_git_blob_sha1')==oid
 (OUT/(oid+'.json')).write_text(json.dumps(proof)+'\n')
 if not proof['verified']:print('FAILED',path,flush=True)
 return proof
result=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
 for i,r in enumerate(pool.map(get,rows)):
  result.append(r)
  if i%50==0:print(i+1,'/',len(rows),'verified',sum(x['verified'] for x in result),flush=True)
(OUT/'result.json').write_text(json.dumps(dict(base=BASE,files=result,all_verified=all(x['verified'] for x in result)),indent=2)+'\n')
assert all(x['verified'] for x in result)
