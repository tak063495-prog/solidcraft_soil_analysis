"""Append only selected typed B-02 statuses from immutable committed bytes."""
import argparse,hashlib,json,pathlib,re,subprocess
p=argparse.ArgumentParser();p.add_argument('--commit',required=True);p.add_argument('names',nargs='+');a=p.parse_args()
assert re.fullmatch('[0-9a-f]{40}',a.commit)
root=pathlib.Path.cwd();folder='interim/2026-12-24-b02-offset-composite-cause-diagnosis'
inventory=root/'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml';auditor=root/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'
before=inventory.read_bytes();tail=''
for name in a.names:
 path=folder+'/'+name+'.json';raw=(root/path).read_bytes();obj=json.loads(raw)
 assert subprocess.check_output(['git','show',a.commit+':'+path])==raw
 sid='b02cause_20261224_'+name;assert ('id = '+json.dumps(sid)) not in before.decode()
 blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
 tail+='\n[[sources]]\nid = '+json.dumps(sid)+'\npath = '+json.dumps(path)+'\ncommit = '+json.dumps(a.commit)+'\nkind = "json"\nsha256 = '+json.dumps(hashlib.sha256(raw).hexdigest())+'\ngit_blob_sha1 = '+json.dumps(blob)+'\n'
 keys=['passed','b02_complete','diagnostic_only']
 if name in ('acceptance','measurements'):keys+=['unsupported_combined_equilibrium']
 for key in keys:
  value=obj[key];assert value is None or type(value) is bool
  tail+='[[sources.status_keys]]\npointer = '+json.dumps('/'+key)+'\nvalue_json = '+json.dumps(json.dumps(value))+'\n'
with inventory.open('ab') as f:f.write(tail.encode())
raw=inventory.read_bytes();assert raw.startswith(before)
s=auditor.read_text(encoding='utf-8');s,count=re.subn(r"SUPPLEMENTAL_SHA256 = '[0-9a-f]{64}'", "SUPPLEMENTAL_SHA256 = '"+hashlib.sha256(raw).hexdigest()+"'",s);assert count==1
with auditor.open('w',encoding='utf-8',newline='\n') as f:f.write(s)
print(json.dumps(dict(added=a.names,commit=a.commit,manifest_sha256=hashlib.sha256(raw).hexdigest())))
