import hashlib,json,pathlib,subprocess,time
root=pathlib.Path.cwd();base='9a5e9d182299b6fce64cc18ea61a65be4702c4e1';start=time.monotonic()
allowed={'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml','interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'}
rows={};differences=[]
for entry in subprocess.check_output(['git','ls-tree','-rz',base]).split(b'\0'):
 if not entry:continue
 meta,name=entry.split(b'\t',1);mode,kind,blob=meta.decode().split();name=name.decode();path=root/name
 if kind!='blob':continue
 h=hashlib.sha1();digest=hashlib.sha256();size=path.stat().st_size;h.update(b'blob '+str(size).encode()+b'\0')
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block);digest.update(block)
 same=h.hexdigest()==blob
 if not same:differences.append(name)
 rows[name]=dict(bytes=size,baseline_git_blob_sha1=blob,current_git_blob_sha1=h.hexdigest(),current_sha256=digest.hexdigest(),identical_to_base=same)
assert not set(differences)-allowed,differences
record=dict(base=base,files=len(rows),unchanged=sum(r['identical_to_base'] for r in rows.values()),changed=sorted(differences),allowed_changes=sorted(allowed),all_other_bytes_identical=True,elapsed_seconds=time.monotonic()-start,sha256_and_git_blobs=rows)
target=root/'work/a03-preserved-files-final.json';assert not target.exists();target.write_text(json.dumps(record,indent=2)+'\n');print({k:v for k,v in record.items() if k!='sha256_and_git_blobs'})
