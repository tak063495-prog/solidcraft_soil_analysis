import hashlib,json,pathlib,re,subprocess,time
root=pathlib.Path('/workspace/solidcraft_soil_analysis');base='1c921b1df0a3da6784476a81c68e6353e39b3678';start=time.monotonic()
allowed={'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml','interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'}
rows={};changed=[]
for entry in subprocess.check_output(['git','ls-tree','-rz',base],cwd=root).split(b'\0'):
 if not entry:continue
 meta,name=entry.split(b'\t',1);mode,kind,blob=meta.decode().split();name=name.decode();p=root/name
 if kind!='blob':continue
 n=p.stat().st_size;h=hashlib.sha1(f'blob {n}\0'.encode());s=hashlib.sha256()
 with p.open('rb') as f:
  while raw:=f.read(1048576):h.update(raw);s.update(raw)
 same=h.hexdigest()==blob
 if not same:changed.append(name)
 rows[name]=dict(size=n,baseline_git_blob_sha1=blob,current_git_blob_sha1=h.hexdigest(),sha256=s.hexdigest(),identical_to_base=same)
assert set(changed)<=allowed,changed
oldinv=subprocess.check_output(['git','show',base+':interim/2026-12-10-ci-p601-inventory-fix/inventory.toml'],cwd=root)
assert (root/'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml').read_bytes().startswith(oldinv)
oldaud=subprocess.check_output(['git','show',base+':interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'],cwd=root)
newaud=(root/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py').read_bytes()
assert re.sub(rb"SUPPLEMENTAL_SHA256 = '[0-9a-f]{64}'",b"SUPPLEMENTAL_SHA256 = 'LEDGER_SHA'",oldaud)==re.sub(rb"SUPPLEMENTAL_SHA256 = '[0-9a-f]{64}'",b"SUPPLEMENTAL_SHA256 = 'LEDGER_SHA'",newaud)
r=dict(base=base,files=len(rows),unchanged=sum(v['identical_to_base'] for v in rows.values()),changed=sorted(changed),all_other_bytes_identical=True,old_ledger_exact_prefix=True,auditor_only_supplemental_sha_constant_changed=True,elapsed_seconds=time.monotonic()-start,files_sha256_and_blob=rows)
out=root/'work/a04cause-preserved-files.json';out.write_text(json.dumps(r,indent=2)+'\n');print({k:v for k,v in r.items() if k!='files_sha256_and_blob'})
