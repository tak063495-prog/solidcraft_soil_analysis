"""Reproduce this bounded diagnosis from the specified commit and sealed evidence.

No production files, checkpoint seals, constitutive equations or tolerances change.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='06bff71dcc54f7060932aff57194f04ffe0593f8'
PARENT='6429e8366ad7b02ab546cc9a9f7ff321b693962b'
sys.path.insert(0,str(ROOT/'tools'))
from verify_high_order import prepare
import probe_nc_midpoint
import analyze
import audit
import path_diagnostic

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)
def prepare_source(out):
    out.mkdir(parents=True,exist_ok=False)
    assert git('rev-parse',BASE+'^').decode().strip()==PARENT
    names=git('ls-tree','-r','--name-only',BASE,'src').decode().splitlines()+['Cargo.toml','Cargo.lock','rust-toolchain.toml']
    for n in names:
        p=out/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(git('show',BASE+':'+n))
    for source,addition in [('mcc.rs','mcc-diagnostic-addition.rs'),('mcc_midpoint.rs','midpoint-diagnostic-addition.rs')]:
        p=out/'src'/source
        p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n')+(HERE/addition).read_bytes())
    (out/'examples').mkdir();shutil.copyfile(HERE/'material_diagnostic.rs',out/'examples/material_diagnostic.rs')
    with (out.parent/(out.name+'-format.log')).open('x') as f:
        subprocess.run(['cargo','fmt','--manifest-path',str(out/'Cargo.toml')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
    return {str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}

def extract_prior(out):
    source=ROOT/'interim/2026-11-09-nc-strict-removal'
    manifest=read(source/'archive-sha256.json')
    if sha(source/'evidence.zip')!=manifest['archive_sha256']:raise ValueError('prior archive SHA256 mismatch')
    with zipfile.ZipFile(source/'evidence.zip') as z:
        if set(z.namelist())!=set(manifest['entries']):raise ValueError('prior archive entry inventory mismatch')
        for n,digest in manifest['entries'].items():
            p=Path(n)
            if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe evidence path')
            b=z.read(n)
            if hashlib.sha256(b).hexdigest()!=digest:raise ValueError('prior entry SHA256 mismatch')
            d=out/p;d.parent.mkdir(parents=True,exist_ok=True);d.write_bytes(b)
    return dict(archive_sha256=manifest['archive_sha256'],verified_entries=len(manifest['entries']))

def prepare_data(out):
    inp=out/'prior/fresh-head-midpoint/input'
    mesh=read(inp/'mesh.json');ground=read(inp/'ground.json');case=read(inp/'stages.json')
    if (mesh,ground,case)!=tuple(probe_nc_midpoint.inputs('midpoint')):raise ValueError('not the specified 5 kPa case')
    q,offsets=prepare(mesh)
    parameters={m['id']:{**m['mcc'],'integration':case['settings']['mcc_integration']} for m in ground['materials']}
    requests=[];records=[]
    for name in ['last-accepted','terminal-rejected']:
        path=out/'prior/experiment/prior/evidence/diagnostic'/name/'trace.jsonl'
        events=[json.loads(l) for l in path.read_text().splitlines()]
        for e in events:
            if e['kind']!='iterate':continue
            label=f'{name}:{e["iteration"]}'
            if max(e['partitions'])>100:raise ValueError('bounded material diagnostic exceeded')
            records.append(dict(label=label,event=e,tangent=next((t for t in events if t['kind']=='tangent' and t['iteration']==e['iteration']),None)))
            for ei,eq in enumerate(q):
                if not e['active'][ei]:continue
                for qi in range(len(eq['b'])):
                    pi=int(offsets[ei]+qi)
                    requests.append(dict(label=label,point=pi,parameters=parameters[mesh['material_ids'][ei]],
                        state=e['committed']['material_states'][pi],strain=e['strain_increment'][pi],
                        partitions=e['partitions'][pi],widths=[1e-7,1e-9,1e-10,1e-11,1e-12]))
    write(out/'material-requests.json',requests);write(out/'field-records.json',records)

def capture(out):
    out.mkdir(parents=True,exist_ok=False)
    write(out/'prior-binding.json',extract_prior(out/'prior'))
    hashes=prepare_source(out/'isolated');write(out/'final-source-sha256.json',hashes)
    target=out/'target';records=[]
    def invoke(args,name,timeout=600):
        args=[str(a) for a in args];t=time.monotonic()
        with (out/(name+'.log')).open('x') as f:
            r=subprocess.run(args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
        records.append(dict(command=args,exit_code=r.returncode,elapsed_seconds=time.monotonic()-t))
        write(out/'capture-processes.json',records);r.check_returncode()
    invoke(['cargo','build','--locked','--release','--manifest-path',out/'isolated/Cargo.toml',
        '--target-dir',target,'--example','material_diagnostic'],'build')
    binary=target/('release/examples/material_diagnostic'+('.exe' if os.name=='nt' else ''))
    write(out/'final-binary-sha256.json',{'material_diagnostic':sha(binary)})
    prepare_data(out)
    invoke([binary,out/'material-requests.json',out/'material-results-final.json'],'material-final',180)
    analyze.measure(out)
    path_diagnostic.run(out,out/'path-diagnostic-coarse',2e-5,64)
    path_diagnostic.run(out,out/'path-diagnostic-fine',1e-5,64)
    result=audit.measure(out,out/'audit')
    after={str(p.relative_to(out/'isolated')):sha(p) for p in sorted((out/'isolated').rglob('*')) if p.is_file()}
    if after!=hashes:raise ValueError('isolated sources changed during diagnosis')
    write(out/'source-after.json',dict(all_exact=True,files=len(hashes)))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source-only',action='store_true');a=p.parse_args()
    if a.source_only:print(json.dumps(prepare_source(a.output.resolve()),indent=2))
    else:print(json.dumps(capture(a.output.resolve())['turning'],indent=2))
