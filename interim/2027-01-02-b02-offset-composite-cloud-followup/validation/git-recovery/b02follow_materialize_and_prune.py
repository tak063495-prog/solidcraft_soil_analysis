from pathlib import Path
import subprocess,os,json,hashlib,shutil
root=Path.cwd();out=root/'work/b02follow-preparation';base='23c7b0ea48906723620b47eefb50e7481074d875'
env=dict(os.environ,GIT_NO_LAZY_FETCH='1');paths={};skip='interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part'
for line in subprocess.check_output(['git','ls-tree','-r',base],text=True,env=env).splitlines():
 meta,path=line.split('\t');mode,kind,oid=meta.split()
 if kind=='blob' and not path.startswith(skip):paths[path]=(oid,mode)
ids=sorted(set(v[0] for v in paths.values()))
r=subprocess.check_output(['git','cat-file','--batch-check=%(objectname) %(objecttype)'],input='\n'.join(ids)+'\n',text=True,env=env)
missing={l.split()[0] for l in r.splitlines() if l.endswith(' missing')};rows=[];written=[]
limit=6*1024**3 if '--initial' in __import__('sys').argv else 100*1024**3
for path,(oid,mode) in sorted(paths.items()):
 p=root/path
 if p.exists() or oid in missing:continue
 if sum(v['bytes'] for v in rows)>limit or shutil.disk_usage(root).free<12*1024**3:break
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('xb') as stream:subprocess.run(['git','cat-file','blob',oid],env=env,stdout=stream,check=True)
 size=p.stat().st_size;h=hashlib.sha1(('blob '+str(size)+'\0').encode());s=hashlib.sha256()
 with p.open('rb') as f:
  for raw in iter(lambda:f.read(1024*1024),b''):h.update(raw);s.update(raw)
 if h.hexdigest()!=oid:raise ValueError(path)
 if mode=='100755':p.chmod(0o755)
 rows.append(dict(path=path,git_blob_sha1=oid,sha256=s.hexdigest(),bytes=size));written.append(path)
 if len(rows)%100==0:print('materialized',len(rows),'free',shutil.disk_usage(root).free,flush=True)
label='initial' if '--initial' in __import__('sys').argv else str(len(list(out.glob('materialization-*.json'))))
(out/('materialization-'+label+'.json')).write_text(json.dumps(dict(base=base,files=rows,free_bytes=shutil.disk_usage(root).free),indent=2)+'\n')
print('complete',len(rows),flush=True)
if '--initial' in __import__('sys').argv:raise SystemExit(0)
verified={}
for report in out.glob('materialization-*.json'):
 for row in json.loads(report.read_text())['files']:verified.setdefault(row['git_blob_sha1'],[]).append(row)
removed=[]
for idx in sorted(Path('.git/objects/pack').glob('pack-*.idx')):
 pack=idx.with_suffix('.pack');promisor=idx.with_suffix('.promisor')
 if not pack.exists() or not promisor.exists() or pack.stat().st_size<1024*1024:continue
 with idx.open('rb') as f:lines=subprocess.check_output(['git','show-index'],stdin=f,text=True).splitlines()
 if len(lines)!=1:continue
 oid=lines[0].split()[1]
 if oid not in verified:continue
 # Every deleted cache blob has at least one exact working copy checked anew.
 for row in verified[oid]:
  p=Path(row['path']);h=hashlib.sha1(('blob '+str(p.stat().st_size)+'\0').encode());s=hashlib.sha256()
  with p.open('rb') as f:
   for raw in iter(lambda:f.read(1024*1024),b''):h.update(raw);s.update(raw)
  assert h.hexdigest()==oid and s.hexdigest()==row['sha256']
 digest=hashlib.file_digest(pack.open('rb'),'sha256').hexdigest()
 removed.append(dict(pack=str(pack),bytes=pack.stat().st_size,pack_sha256=digest,oid=oid,verified_paths=[v['path'] for v in verified[oid]]))
 pack.unlink();idx.unlink();promisor.unlink();idx.with_suffix('.rev').unlink(missing_ok=True)
(out/('cache-pruning-'+label+'.json')).write_text(json.dumps(dict(working_files_or_commits_deleted=False,only_verified_refetchable_single_blob_cache=True,removed=removed,free_bytes=shutil.disk_usage('.').free),indent=2)+'\n')
print('pruned',len(removed),flush=True)
