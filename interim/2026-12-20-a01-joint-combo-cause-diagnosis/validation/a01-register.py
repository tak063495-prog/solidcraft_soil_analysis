import hashlib,json,re,subprocess,tomllib
from pathlib import Path
root=Path('/workspace/work/a01-source');commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
manifest=root/'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml';before=manifest.read_bytes();data=tomllib.loads(before.decode());used={s['path'] for s in data['sources']}
number=max(int(s['id'][1:]) for s in data['sources'] if re.fullmatch(r'N\d+',s['id']))+1
pieces=[]
for filename in ('acceptance.json','measurements.json'):
 path='interim/2026-12-20-a01-joint-combo-cause-diagnosis/'+filename
 assert path not in used
 raw=(root/path).read_bytes();value=json.loads(raw)
 blob=subprocess.check_output(['git','rev-parse',commit+':'+path],cwd=root,text=True).strip()
 assert subprocess.check_output(['git','show',commit+':'+path],cwd=root)==raw
 keys=[]
 def walk(value,pointer=''):
  if value is None or type(value) is bool:keys.append((pointer,value))
  elif isinstance(value,dict):
   for k,v in sorted(value.items()):walk(v,pointer+'/'+k.replace('~','~0').replace('/','~1'))
  elif isinstance(value,list):
   for i,v in enumerate(value):walk(v,pointer+'/'+str(i))
 walk(value);assert keys
 text=f'\n[[sources]]\nid = "N{number:03}"\npath = "{path}"\ncommit = "{commit}"\nkind = "json"\nsha256 = "{hashlib.sha256(raw).hexdigest()}"\ngit_blob_sha1 = "{blob}"\n'
 for pointer,value in keys:text+='[[sources.status_keys]]\npointer = '+json.dumps(pointer)+'\nvalue_json = '+json.dumps(json.dumps(value))+'\n'
 pieces.append(text);number+=1
manifest.write_bytes(before+''.join(pieces).encode());assert manifest.read_bytes().startswith(before)
digest=hashlib.sha256(manifest.read_bytes()).hexdigest();audit=root/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py';s=audit.read_text();s,n=re.subn(r"SUPPLEMENTAL_SHA256 = '[a-f0-9]{64}'","SUPPLEMENTAL_SHA256 = '"+digest+"'",s);assert n==1;audit.write_text(s)
print(commit,digest)
