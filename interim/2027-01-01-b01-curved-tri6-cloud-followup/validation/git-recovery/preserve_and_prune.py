"""Verify all base files, retain permitted ledger edits, then prune archive cache."""
from pathlib import Path
import subprocess,hashlib,json,os,shutil
root=Path.cwd();out=root/'work/b01follow-git-recovery';base='c6b1176b1147b8f242ef979f6c4a6a70e6829579';env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
allowed={'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py','interim/2026-12-10-ci-p601-inventory-fix/inventory.toml'}
paths=[]
for line in subprocess.check_output(['git','ls-tree','-r',base],text=True,env=env).splitlines():
 meta,path=line.split('\t');mode,kind,oid=meta.split()
 if kind=='blob':paths.append((path,mode,oid))
verified=[];excluded=[];changed=[];archives={};total=0
for path,mode,oid in paths:
 p=root/path
 if path.startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part'):
  excluded.append(dict(path=path,git_blob_sha1=oid,existing_ci_exclusion=True));continue
 if not p.is_file():raise ValueError('missing base file '+path)
 size=p.stat().st_size;h=hashlib.sha1(('blob '+str(size)+'\0').encode());s=hashlib.sha256()
 with p.open('rb') as f:
  for raw in iter(lambda:f.read(1024*1024),b''):h.update(raw);s.update(raw)
 row=dict(path=path,git_blob_sha1=oid,actual_blob_sha1=h.hexdigest(),sha256=s.hexdigest(),bytes=size,
          executable_bit_matches=(bool(p.stat().st_mode&0o111)==(mode=='100755')))
 if not row['executable_bit_matches']:raise ValueError('mode '+path)
 if h.hexdigest()!=oid:
  if path not in allowed:raise ValueError('unapproved change '+path)
  old=subprocess.check_output(['git','show',base+':'+path],env=env);current=p.read_bytes()
  if path.endswith('inventory.toml'):
   assert current.startswith(old);row['old_bytes_exact_prefix']=True
  else:
   import re
   replace=lambda raw:re.sub(rb"SUPPLEMENTAL_SHA256 = '[0-9a-f]{64}'",b"SUPPLEMENTAL_SHA256 = 'PIN'",raw)
   assert replace(old)==replace(current);row['only_inventory_pin_changed']=True
  changed.append(row)
 else:
  verified.append(row);total+=size
  if any(suffix in path for suffix in ('.zip','.tar','.gz','.npz','.part')):archives.setdefault(oid,[]).append(row)
 if (len(verified)+len(changed))%1000==0:print('verified',len(verified),'bytes',total,flush=True)
assert {r['path'] for r in changed}==allowed
preservation=dict(base=base,base_paths=len(paths),unchanged_files=len(verified),unchanged_bytes=total,
 unchanged=verified,permitted_ledger_edits=changed,unmaterialized_existing_exclusions=excluded,old_numeric_evidence_and_failures_unchanged=True,
 production_rust_and_nc_unchanged=True,forbidden_changes=[])
(out/'base-preservation.json').write_text(json.dumps(preservation,indent=2)+'\n')
removed=[];before=shutil.disk_usage(root).free
for idx in sorted((root/'.git/objects/pack').glob('pack-*.idx')):
 pack=idx.with_suffix('.pack');promisor=idx.with_suffix('.promisor')
 if not pack.exists() or not promisor.exists() or pack.stat().st_size<1024*1024:continue
 with idx.open('rb') as f:lines=subprocess.check_output(['git','show-index'],stdin=f,text=True).splitlines()
 if len(lines)!=1:continue
 oid=lines[0].split()[1]
 if oid not in archives:continue
 with pack.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 removed.append(dict(pack=pack.name,bytes=pack.stat().st_size,pack_sha256=digest,git_blob_sha1=oid,verified_paths=[v['path'] for v in archives[oid]]))
 pack.unlink();idx.unlink();promisor.unlink();idx.with_suffix('.rev').unlink(missing_ok=True)
(out/'cache-pruning-final.json').write_text(json.dumps(dict(only_refetchable_single_archive_promisor_cache=True,working_files_or_commits_deleted=False,
 removed=removed,freed_bytes=shutil.disk_usage(root).free-before,free_bytes=shutil.disk_usage(root).free),indent=2)+'\n')
print('all base bytes preserved; free',shutil.disk_usage(root).free,flush=True)
