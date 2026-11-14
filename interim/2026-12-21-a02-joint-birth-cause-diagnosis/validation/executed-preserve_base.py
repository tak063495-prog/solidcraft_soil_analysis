import pathlib,subprocess,json,time,hashlib,re
root=pathlib.Path('/workspace/solidcraft_soil_analysis');out=root/'interim/2026-12-21-a02-joint-birth-cause-diagnosis';base='7b254b1cdbe9c288989631944cd698aa861f83ac'
while not (root/'work/a02-materialization.json').exists():time.sleep(1)
rows=[];changed=[];missing=[];by={}
allowed={'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py','interim/2026-12-10-ci-p601-inventory-fix/inventory.toml'}
for line in subprocess.check_output(['git','ls-tree','-r',base],cwd=root,text=True).splitlines():
 meta,name=line.split('\t');mode,kind,oid=meta.split()
 if kind!='blob':continue
 path=root/name
 if not path.is_file():missing.append(name);continue
 size=path.stat().st_size;blob=hashlib.sha1(b'blob '+str(size).encode()+b'\0');sha=hashlib.sha256()
 with path.open('rb') as f:
  for raw in iter(lambda:f.read(1024*1024),b''):blob.update(raw);sha.update(raw)
 same=blob.hexdigest()==oid
 rows.append(dict(path=name,base_git_blob_sha1=oid,working_sha256=sha.hexdigest(),unchanged=same))
 group=name.split('/')[0];g=by.setdefault(group,dict(files=0,unchanged=0));g['files']+=1;g['unchanged']+=same
 if not same:changed.append(name)
assert not missing and set(changed)==allowed,(missing,changed)
checks={}
for name in allowed:
 raw=subprocess.check_output(['git','show',base+':'+name],cwd=root);now=(root/name).read_bytes()
 if name.endswith('inventory.toml'):checks['old_inventory_is_exact_prefix']=now.startswith(raw)
 else:
  pattern=rb"SUPPLEMENTAL_SHA256 = '[a-f0-9]{64}'"
  checks['audit_changed_only_inventory_digest']=re.sub(pattern,b'SUPPLEMENTAL_SHA256 = PINNED',raw)==re.sub(pattern,b'SUPPLEMENTAL_SHA256 = PINNED',now)
assert all(checks.values())
result=dict(base_commit=base,checked_at_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
 base_tracked_blobs=len(rows),unchanged=len(rows)-len(changed),changed_existing_paths=changed,missing_paths=missing,
 by_root=by,all_original_paths_present=True,all_protected_numeric_evidence_unchanged=True,
 production_rust_unchanged=True,nc_code_and_configuration_unchanged=True,old_hash_snapshot_and_failures_unchanged=True,registration_only_checks=checks)
(out/'base-preservation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
(root/'work/a02-validation/base-byte-manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
print(result)
