import concurrent.futures,hashlib,json,subprocess,sys,time,urllib.parse
from pathlib import Path
root=Path.cwd();base=root/'work/extra-artifacts';base.mkdir(exist_ok=True)
tree={}
for line in subprocess.check_output(['git','ls-tree','-r','HEAD'],text=True).splitlines():
 meta,name=line.split('\t');mode,kind,oid=meta.split()
 if kind=='blob':tree[name]=oid
paths=sys.argv[1:]
def get(name):
 assert name in tree,name
 target=root/name;oid=tree[name];temp=base/(oid+'.download')
 if target.exists():raw=target.read_bytes();source='existing'
 else:
  target.parent.mkdir(parents=True,exist_ok=True)
  url='https://raw.githubusercontent.com/tak063495-prog/solidcraft_soil_analysis/3b1133d3be9905dc24f62c27b5b1eb456a83f2e2/'+urllib.parse.quote(name)
  with (base/(oid+'.log')).open('w') as f:
   r=subprocess.run(['curl','-fSL','--retry','2','--max-time','240',url,'-o',str(temp)],stdout=f,stderr=subprocess.STDOUT)
  assert r.returncode==0,(name,r.returncode)
  raw=temp.read_bytes();source='original blob from fixed-commit URL'
 assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==oid,name
 if not target.exists():temp.rename(target)
 return dict(path=name,git_blob=oid,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),source=source)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:values=list(pool.map(get,paths))
p=base/('records-'+str(time.time_ns())+'.json');p.write_text(json.dumps(values,indent=2)+'\n');print('verified',len(values),p)
