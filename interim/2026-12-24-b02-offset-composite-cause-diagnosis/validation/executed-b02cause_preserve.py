import hashlib,json,pathlib,re,subprocess,time
root=pathlib.Path(__file__).resolve().parents[1];base='33b185b060c26cf2d0957271391b7bec86e491da';start=time.monotonic();rows=[]
for line in subprocess.check_output(['git','ls-tree','-r',base],cwd=root,text=True).splitlines():
 info,name=line.split('\t',1);mode,kind,expected=info.split()
 if kind!='blob':continue
 p=root/name
 if not p.is_file():rows.append(dict(path=name,missing=True));continue
 h=hashlib.sha1(b'blob '+str(p.stat().st_size).encode()+b'\0');s=hashlib.sha256()
 with p.open('rb') as f:
  for raw in iter(lambda:f.read(1024*1024),b''):h.update(raw);s.update(raw)
 rows.append(dict(path=name,size=p.stat().st_size,expected_git_blob=expected,actual_git_blob=h.hexdigest(),sha256=s.hexdigest(),same=h.hexdigest()==expected,executable_mode_matches=bool(p.stat().st_mode&0o111)==(mode=='100755')))
allowed={'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml','interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'}
missing=[r['path'] for r in rows if r.get('missing')];changed=[r['path'] for r in rows if not r.get('same',False)]
assert not missing and set(changed)==allowed
inv='interim/2026-12-10-ci-p601-inventory-fix/inventory.toml';old_inv=subprocess.check_output(['git','show',base+':'+inv],cwd=root);assert (root/inv).read_bytes().startswith(old_inv)
audit='interim/2026-12-06-p6-01-status-reconciliation/audit_status.py';old_audit=subprocess.check_output(['git','show',base+':'+audit],cwd=root)
mask=lambda raw:re.sub(rb"SUPPLEMENTAL_SHA256 = '[a-f0-9]{64}'",b"SUPPLEMENTAL_SHA256 = 'MASKED_DIGEST_ONLY'",raw)
assert mask(old_audit)==mask((root/audit).read_bytes())
assert all(r['executable_mode_matches'] for r in rows)
result=dict(schema='b02cause-preservation/1',diagnostic_only=True,b02_complete=False,integrity_ok=True,base_commit=base,baseline_files=len(rows),unchanged_files=len(rows)-len(changed),changed_existing_paths=changed,missing=missing,executable_modes_match=True,inventory_original_bytes_retained_as_prefix=True,auditor_digest_constant_only=True,old_numeric_values_or_seals_changed=False,elapsed_seconds=time.monotonic()-start,files=rows)
(root/'interim/2026-12-24-b02-offset-composite-cause-diagnosis/preservation.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k!='files'})
