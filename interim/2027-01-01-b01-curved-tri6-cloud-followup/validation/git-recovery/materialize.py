"""Materialize available fixed blobs and verify original bytes, bounded by free disk."""
from pathlib import Path
import subprocess,os,json,hashlib,shutil
root=Path.cwd();out=root/'work/b01follow-git-recovery';base='c6b1176b1147b8f242ef979f6c4a6a70e6829579'
env=dict(os.environ,GIT_NO_LAZY_FETCH='1');paths={}
for line in subprocess.check_output(['git','ls-tree','-r',base],text=True,env=env).splitlines():
 meta,path=line.split('\t');mode,kind,oid=meta.split()
 if kind=='blob' and not path.startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part'):
  paths[path]=(oid,mode)
ids=sorted(set(v[0] for v in paths.values()))
r=subprocess.check_output(['git','cat-file','--batch-check=%(objectname) %(objecttype)'],input='\n'.join(ids)+'\n',text=True,env=env)
missing={l.split()[0] for l in r.splitlines() if l.endswith(' missing')};rows=[];written=[]
for path,(oid,mode) in sorted(paths.items()):
 p=root/path
 if p.exists() or oid in missing:continue
 if shutil.disk_usage(root).free<25*1024**3:break
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('xb') as stream:subprocess.run(['git','cat-file','blob',oid],env=env,stdout=stream,check=True)
 size=p.stat().st_size;h=hashlib.sha1(('blob '+str(size)+'\0').encode());s=hashlib.sha256()
 with p.open('rb') as f:
  for raw in iter(lambda:f.read(1024*1024),b''):h.update(raw);s.update(raw)
 if h.hexdigest()!=oid:raise ValueError(path)
 if mode=='100755':p.chmod(0o755)
 rows.append(dict(path=path,git_blob_sha1=oid,sha256=s.hexdigest(),bytes=size,original_blob_matches=True));written.append(path)
 if len(rows)%100==0:print('materialized',len(rows),'free',shutil.disk_usage(root).free,flush=True)
if written:
 subprocess.run(['git','update-index','--no-skip-worktree','--stdin'],input='\n'.join(written)+'\n',text=True,check=True)
 subprocess.run(['git','update-index','--refresh'],stdout=subprocess.DEVNULL,check=False)
(out/'materialization.json').write_text(json.dumps(dict(base=base,files=rows,initial_missing_objects=len(missing),free_bytes=shutil.disk_usage(root).free),indent=2)+'\n')
print('complete',len(rows),flush=True)
