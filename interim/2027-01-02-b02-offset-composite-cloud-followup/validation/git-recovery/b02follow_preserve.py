import subprocess,hashlib,json,os,re
from pathlib import Path
base='23c7b0ea48906723620b47eefb50e7481074d875';root=Path.cwd();D=root/'interim/2027-01-02-b02-offset-composite-cloud-followup'
rows=[];errors=[];excluded=[];special=[]
inv='interim/2026-12-10-ci-p601-inventory-fix/inventory.toml';audit='interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'
env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
for line in subprocess.check_output(['git','ls-tree','-r',base],text=True,env=env).splitlines():
 meta,name=line.split('\t');mode,kind,oid=meta.split()
 if kind!='blob':continue
 if name.startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part'):
  excluded.append(dict(path=name,oid=oid,git_reference_preserved=True,existing_ci_exclusion=True));continue
 p=root/name
 if not p.is_file():errors.append(dict(path=name,error='missing'));continue
 if name in (inv,audit):
  old=subprocess.check_output(['git','show',base+':'+name],env=env);new=p.read_bytes()
  if name==inv:ok=new.startswith(old)
  else:
   strip=lambda t:re.sub(rb"SUPPLEMENTAL_SHA256 = '[0-9a-f]{64}'",b"SUPPLEMENTAL_SHA256 = 'PIN'",t)
   ok=strip(old)==strip(new)
  special.append(dict(path=name,allowed_change_only=ok,old_sha256=hashlib.sha256(old).hexdigest(),new_sha256=hashlib.sha256(new).hexdigest()))
  if not ok:errors.append(dict(path=name,error='unallowed special edit'))
  continue
 size=p.stat().st_size;h=hashlib.sha1(('blob '+str(size)+'\0').encode());s=hashlib.sha256()
 with p.open('rb') as f:
  for raw in iter(lambda:f.read(1024*1024),b''):h.update(raw);s.update(raw)
 ok=h.hexdigest()==oid and ((p.stat().st_mode&0o111)!=0)==(mode=='100755')
 rows.append(dict(path=name,git_blob_sha1=oid,sha256=s.hexdigest(),bytes=size,exact=ok))
 if not ok:errors.append(dict(path=name,error='bytes/mode differ'))
 if len(rows)%1500==0:print('checked',len(rows),flush=True)
result=dict(base=base,exact_old_files=len(rows),bytes=sum(r['bytes'] for r in rows),special=special,existing_ci_exclusions=excluded,errors=errors,all_preserved=not errors)
(D/'base-tree-binding.json').write_text(json.dumps(rows,indent=2)+'\n')
Path(__import__('sys').argv[1]).write_text(json.dumps(result,indent=2)+'\n')
print('preserved',len(rows),'errors',len(errors),flush=True)
if errors:raise SystemExit(1)
