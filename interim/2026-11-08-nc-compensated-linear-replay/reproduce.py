"""Bounded reproduction only: saved systems, controls, and one minimum increment."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import zipfile
from prepare import prepare, ROOT, HERE
from audit_linear import audit
from measure_trace import measure


def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def write(p, v):
    with p.open('x') as f: json.dump(v, f, indent=2); f.write('\n')

def read(p): return json.loads(p.read_text())


def matrices(path, label):
    return [dict(label=label+':'+str(e['iteration']),matrix_columns=e['matrix_columns'],reduced_residual=e['reduced_residual'])
            for e in map(json.loads,path.read_text().splitlines()) if e['kind']=='tangent']


def stripped(v):
    if isinstance(v,dict):return {k:stripped(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stripped(x) for x in v]
    return v


def capture(out):
    out.mkdir(parents=True,exist_ok=False)
    prior=ROOT/'interim/2026-11-08-nc-midpoint-stop-diagnosis'
    manifest=read(prior/'archive-sha256.json');assert digest(prior/'evidence.zip')==manifest['archive_sha256']
    with zipfile.ZipFile(prior/'evidence.zip') as z:
        for n,h in manifest['entries'].items():
            p=Path(n);assert not p.is_absolute() and '..' not in p.parts
            b=z.read(n);assert hashlib.sha256(b).hexdigest()==h
            dest=out/'evidence'/n;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(b)
    write(out/'prior-binding.json',manifest)
    records=[]
    def run(args,name):
        args=[str(a) for a in args];start=time.monotonic()
        with (out/(name+'.log')).open('x') as f:p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,timeout=300)
        records.append(dict(command=args,exit_code=p.returncode,elapsed_seconds=time.monotonic()-start))
        (out/'processes.json').write_text(json.dumps(records,indent=2)+'\n');p.check_returncode()
    ext='.exe' if os.name=='nt' else ''
    binaries={}
    for name,candidate in [('baseline',False),('candidate',True)]:
        source=out/name;prepare(source,candidate)
        run(['cargo','fmt','--manifest-path',source/'Cargo.toml'],name+'-fmt')
        run(['cargo','build','--locked','--release','--manifest-path',source/'Cargo.toml','--example','nc_stop_probe','--example','linear_probe'],name+'-build')
        binaries[name]=source/('target/release/examples/nc_stop_probe'+ext)
    write(out/'binary-sha256.json',{k:digest(v) for k,v in binaries.items()})
    req=[]
    for n in ('last-accepted','terminal-rejected'):req+=matrices(out/'evidence/diagnostic'/n/'trace.jsonl',n)
    write(out/'matrices.json',req)
    linear=out/('candidate/target/release/examples/linear_probe'+ext)
    run([linear,out/'matrices.json',out/'linear-results.json'],'linear-probe')
    linear_audit=audit(req,read(out/'linear-results.json'));write(out/'linear-audit.json',linear_audit)
    # The nonlinear trial is forbidden until every saved system passes independently.
    if not linear_audit['all_saved_matrices_pass_independent_linear_gate']:
        return dict(passed=False,history_converged=False,nonlinear_trial_performed=False,reason='independent linear gate failed')
    folder=out/'evidence/prior/nc-resumed'
    mesh=folder/'input/mesh.json'
    for n,b,cp,target in [('baseline-trial','baseline',0,'0.9514002799987793'),
                           ('candidate-accepted-control','candidate',1,'0.9514000415802002'),
                           ('candidate-trial','candidate',0,'0.9514002799987793')]:
        run([binaries[b],mesh,folder/'job'/f'recovery-checkpoint-{cp}.json',target,out/n],n)
    controls=dict(last_accepted_state_exact=read(out/'candidate-accepted-control/diagnostic.json')['outcome'].get('state')==read(folder/'job/recovery-checkpoint-0.json')['run']['last_committed'],
                  baseline_trace_exact_except_elapsed=stripped([json.loads(l) for l in (out/'baseline-trial/trace.jsonl').read_text().splitlines()])==stripped([json.loads(l) for l in (out/'evidence/diagnostic/terminal-rejected/trace.jsonl').read_text().splitlines()]))
    write(out/'controls.json',controls);assert all(controls.values())
    req=matrices(out/'candidate-trial/trace.jsonl','candidate-trial');write(out/'candidate-matrices.json',req)
    run([linear,out/'candidate-matrices.json',out/'candidate-linear-results.json'],'candidate-linear-probe')
    write(out/'candidate-linear-audit.json',audit(req,read(out/'candidate-linear-results.json')))
    result=measure(out/'evidence',out/'candidate-trial',out/'measurements')
    write(out/'bounded-result.json',dict(passed=False,history_converged=False,full_removal_reached=False,
        target_trial_converged=read(out/'candidate-trial/diagnostic.json')['outcome']['converged'],
        accepted_new_increments=result['accepted_new_increments']))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=capture(a.output.resolve());print(json.dumps({k:v for k,v in result.items() if k not in ('force_rows','spectra')},indent=2))
