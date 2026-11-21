"""Append-only A-04 P6 registration with fixed commit and typed status proofs."""
import argparse,hashlib,json,pathlib,re,subprocess
p=argparse.ArgumentParser();p.add_argument('--commit',required=True);p.add_argument('names',nargs='+');a=p.parse_args()
root=pathlib.Path.cwd();folder='interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis'
inventory=root/'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml';auditor=root/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'
before=inventory.read_bytes();tail=''
for name in a.names:
 path=folder+'/'+name+'.json';raw=(root/path).read_bytes();obj=json.loads(raw)
 fixed=subprocess.check_output(['git','show',a.commit+':'+path]);assert fixed==raw
 sid='a04cause_20261223_'+name;assert ('id = '+json.dumps(sid)) not in before.decode()
 blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
 tail+='\n[[sources]]\nid = '+json.dumps(sid)+'\npath = '+json.dumps(path)+'\ncommit = '+json.dumps(a.commit)+'\nkind = "json"\nsha256 = '+json.dumps(hashlib.sha256(raw).hexdigest())+'\ngit_blob_sha1 = '+json.dumps(blob)+'\n'
 proofs={k:obj[k] for k in ('passed','a04_complete','diagnostic_only')}
 if name=='measurements':proofs.update({'interface/'+k:v for k,v in obj['interface'].items()})
 if name=='acceptance':proofs.update({k:obj[k] for k in ('interface_head_mapping','physical_time_generation','joint_transient_water_balance','coupled_time_refinement')})
 for key,value in proofs.items():
  assert value is None or type(value) is bool
  tail+='[[sources.status_keys]]\npointer = '+json.dumps('/'+key)+'\nvalue_json = '+json.dumps(json.dumps(value))+'\n'
with inventory.open('ab') as f:f.write(tail.encode())
raw=inventory.read_bytes();assert raw.startswith(before)
s=auditor.read_text(encoding='utf-8');s,count=re.subn(r"SUPPLEMENTAL_SHA256 = '[0-9a-f]{64}'", "SUPPLEMENTAL_SHA256 = '"+hashlib.sha256(raw).hexdigest()+"'",s);assert count==1
with auditor.open('w',encoding='utf-8',newline='\n') as f:f.write(s)
print(json.dumps(dict(added=a.names,commit=a.commit,manifest_sha256=hashlib.sha256(raw).hexdigest())))
