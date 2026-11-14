"""Append-only P6 source registration after concrete evidence commits exist."""
import argparse,hashlib,json,pathlib,re,subprocess
p=argparse.ArgumentParser();p.add_argument('--commit',required=True);p.add_argument('names',nargs='+');a=p.parse_args()
root=pathlib.Path.cwd();folder='interim/2026-12-22-a03-wet-soil-birth-cause-diagnosis'
inventory=root/'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml';auditor=root/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'
with inventory.open('rb') as f:before=f.read()
tail=''
for name in a.names:
 path=folder+'/'+name+'.json';raw=(root/path).read_bytes();obj=json.loads(raw)
 fixed=subprocess.check_output(['git','show',a.commit+':'+path]);assert fixed==raw
 sid='a03cause_20261222_'+name;assert ('id = '+json.dumps(sid)) not in before.decode()
 blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
 tail+='\n[[sources]]\nid = '+json.dumps(sid)+'\npath = '+json.dumps(path)+'\ncommit = '+json.dumps(a.commit)+'\nkind = "json"\nsha256 = '+json.dumps(hashlib.sha256(raw).hexdigest())+'\ngit_blob_sha1 = '+json.dumps(blob)+'\n'
 for key in ('passed','a03_complete','diagnostic_only'):
  if key not in obj:continue
  tail+='[[sources.status_keys]]\npointer = '+json.dumps('/'+key)+'\nvalue_json = '+json.dumps(json.dumps(obj[key]))+'\n'
 if name=='measurements':
  for key in ('flux','storage','birth_water_transfer','balance_acceptance'):
   assert obj['water_accounting'][key] is None
   tail+='[[sources.status_keys]]\npointer = '+json.dumps('/water_accounting/'+key)+'\nvalue_json = "null"\n'
with inventory.open('ab') as f:f.write(tail.encode())
raw=inventory.read_bytes();assert raw.startswith(before)
s=auditor.read_text(encoding='utf-8');s,count=re.subn(r"SUPPLEMENTAL_SHA256 = '[0-9a-f]{64}'", "SUPPLEMENTAL_SHA256 = '"+hashlib.sha256(raw).hexdigest()+"'",s);assert count==1
with auditor.open('w',encoding='utf-8',newline='\n') as f:f.write(s)
print(json.dumps(dict(added=a.names,commit=a.commit,manifest_sha256=hashlib.sha256(raw).hexdigest())))
