"""Append fixed-commit P6 sources, keeping the complete old inventory prefix."""
import argparse, hashlib, json, os
from pathlib import Path
import re, subprocess

ROOT=Path(__file__).resolve().parents[2]
PREFIX='interim/2026-12-28-a01-joint-combo-cloud-followup/'
INVENTORY=ROOT/'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml'
AUDITOR=ROOT/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'

def main():
    p=argparse.ArgumentParser();p.add_argument('--commit',required=True);p.add_argument('--source',choices=('scientific','validation'),required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if not re.fullmatch('[0-9a-f]{40}',a.commit):raise ValueError('fixed full SHA required')
    env=dict(os.environ,GIT_NO_LAZY_FETCH='1');before=INVENTORY.read_bytes();addition='';proofs=[]
    for name in (('acceptance','measurements') if a.source=='scientific' else ('validation',)):
        path=PREFIX+name+'.json';raw=(ROOT/path).read_bytes()
        fixed=subprocess.check_output(['git','show',a.commit+':'+path],cwd=ROOT,env=env)
        if raw!=fixed:raise ValueError('working source differs from fixed commit')
        value=json.loads(raw);sid='a01follow_20261228_'+name
        if ('id = "'+sid+'"').encode() in before:raise ValueError('duplicate registration')
        digest=hashlib.sha256(raw).hexdigest();blob=hashlib.sha1(('blob '+str(len(raw))+'\0').encode()+raw).hexdigest()
        actual=subprocess.check_output(['git','rev-parse',a.commit+':'+path],cwd=ROOT,env=env,text=True).strip()
        if actual!=blob:raise ValueError('Git blob differs')
        addition+='\n[[sources]]\n'+''.join(k+' = '+json.dumps(v)+'\n' for k,v in
            dict(id=sid,path=path,commit=a.commit,kind='json',sha256=digest,git_blob_sha1=blob).items())
        for key in ('passed','a01_complete','diagnostic_only'):
            if type(value[key]) is not bool:raise ValueError('status must be typed bool')
            addition+='[[sources.status_keys]]\npointer = '+json.dumps('/'+key)+'\nvalue_json = '+json.dumps(json.dumps(value[key]))+'\n'
        proofs.append(dict(path=path,commit=a.commit,sha256=digest,git_blob_sha1=blob))
    after=before+addition.encode('utf-8');INVENTORY.write_bytes(after)
    old_hash=hashlib.sha256(before).hexdigest();new_hash=hashlib.sha256(after).hexdigest()
    source=AUDITOR.read_text(encoding='utf-8');old="SUPPLEMENTAL_SHA256 = '"+old_hash+"'"
    if source.count(old)!=1:raise ValueError('auditor pin mismatch')
    updated=source.replace(old,"SUPPLEMENTAL_SHA256 = '"+new_hash+"'")
    AUDITOR.write_text(updated,encoding='utf-8',newline='\n')
    with a.output.open('x') as f:json.dump(dict(old_inventory_exact_prefix=after.startswith(before),new_sources=proofs,
        inventory_sha256=new_hash,auditor_only_pin_changed=source.replace(old,'PIN')==updated.replace("SUPPLEMENTAL_SHA256 = '"+new_hash+"'",'PIN')),f,indent=2);f.write('\n')

if __name__=='__main__':main()
