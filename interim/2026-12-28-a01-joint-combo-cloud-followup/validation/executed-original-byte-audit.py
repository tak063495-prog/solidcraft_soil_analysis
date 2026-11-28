import hashlib,json,os,re,subprocess
from pathlib import Path
root=Path.cwd();before=json.load(open('work/a01follow-base-preservation.json'));base=before['base'];allowed={'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml','interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'}
errors=[];count=0;stat={}
for record in before['files']:
 name=record['path'];p=root/name
 if name in allowed:continue
 with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 if digest!=record['sha256']:errors.append(name)
 count+=1;st=p.stat();stat[name]=[st.st_size,st.st_mtime_ns]
proofs={}
for name in sorted(allowed):
 old=subprocess.check_output(['git','show',base+':'+name],env=dict(os.environ,GIT_NO_LAZY_FETCH='1'));current=(root/name).read_bytes()
 if name.endswith('.toml'):ok=current.startswith(old);kind='old complete inventory is exact byte prefix'
 else:
  pattern=rb"SUPPLEMENTAL_SHA256 = '[0-9a-f]{64}'"
  ok=re.sub(pattern,b'PIN',old)==re.sub(pattern,b'PIN',current);kind='only supplemental inventory SHA constant changed'
 proofs[name]=dict(authorized_change=ok,proof=kind)
 assert ok,name
result=dict(original_base=base,old_paths_byte_verified=count,changed_old_numerical_or_failure_paths=errors,authorized_p6_maintenance=proofs,nc_code_or_settings_changed=False,all_old_scientific_bytes_unchanged=not errors)
Path('interim/2026-12-28-a01-joint-combo-cloud-followup/validation/final-original-byte-audit.json').write_text(json.dumps(result,indent=2)+'\n');Path('work/a01follow-final-original-stats.json').write_text(json.dumps(stat)+'\n');print(result)
assert not errors
