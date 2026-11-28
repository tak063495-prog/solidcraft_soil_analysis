"""Register immutable scientific/validation records using existing P6 SHA scheme."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
ROOT=Path(__file__).resolve().parents[2]
PREFIX='interim/2027-01-08-pc-failure-diagnosis/'
INVENTORY=ROOT/'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml'
AUDITOR=ROOT/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'


def register(commit,kind,output):
    if not re.fullmatch('[0-9a-f]{40}',commit):raise ValueError('fixed full SHA required')
    names=('acceptance','measurements') if kind=='scientific' else ('validation',)
    env=dict(os.environ,GIT_NO_LAZY_FETCH='1');before=INVENTORY.read_bytes();addition='';proofs=[]
    for name in names:
        path=PREFIX+name+'.json';raw=(ROOT/path).read_bytes()
        fixed=subprocess.check_output(['git','show',commit+':'+path],cwd=ROOT,env=env)
        if raw!=fixed:raise ValueError('source differs from fixed commit')
        value=json.loads(raw);sid='pcfail_20270108_'+name
        if ('id = "'+sid+'"').encode() in before:raise ValueError('duplicate registration')
        digest=hashlib.sha256(raw).hexdigest();blob=hashlib.sha1(('blob '+str(len(raw))+'\0').encode()+raw).hexdigest()
        actual=subprocess.check_output(['git','rev-parse',commit+':'+path],cwd=ROOT,env=env,text=True).strip()
        if actual!=blob:raise ValueError('Git blob differs')
        addition+='\n[[sources]]\n'+''.join(k+' = '+json.dumps(v)+'\n' for k,v in
            dict(id=sid,path=path,commit=commit,kind='json',sha256=digest,git_blob_sha1=blob).items())
        for key in ('passed','diagnostic_only','all_acceptance_complete'):
            if type(value[key]) is not bool:raise ValueError('typed status bool required')
            addition+='[[sources.status_keys]]\npointer = '+json.dumps('/'+key)+'\nvalue_json = '+json.dumps(json.dumps(value[key]))+'\n'
        proofs.append(dict(path=path,commit=commit,sha256=digest,git_blob_sha1=blob))
    after=before+addition.encode();oldhash=hashlib.sha256(before).hexdigest();newhash=hashlib.sha256(after).hexdigest()
    source=AUDITOR.read_text(encoding='utf-8');old="SUPPLEMENTAL_SHA256 = '"+oldhash+"'"
    if source.count(old)!=1:raise ValueError('auditor pin differs')
    updated=source.replace(old,"SUPPLEMENTAL_SHA256 = '"+newhash+"'")
    # Validate both mutations before writing either.
    INVENTORY.write_bytes(after);AUDITOR.write_text(updated,encoding='utf-8',newline='\n')
    with output.open('x',encoding='utf-8') as f:
        json.dump(dict(old_inventory_exact_prefix=after.startswith(before),new_sources=proofs,inventory_sha256=newhash,
                      auditor_only_pin_changed=source.replace(old,'PIN')==updated.replace("SUPPLEMENTAL_SHA256 = '"+newhash+"'",'PIN')),f,indent=2);f.write('\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--commit',required=True);p.add_argument('--kind',choices=('scientific','validation'),required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    register(a.commit,a.kind,a.output)
