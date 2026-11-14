"""Reproduce the bounded offset acceptance and retain curved-shell refusals."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
from plate6_offset_study import study,write,zero_equivalence,ROOT
from plate6_scope_rejections import probe

def record(output,binary=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    binary=Path(binary or ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d'))
    cmd=['cargo','test','--locked','--release','--test','plates6_offset','--','--nocapture']
    kernel=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True)
    (output/'kernel.log').write_text(kernel.stdout+kernel.stderr,encoding='utf-8')
    if kernel.returncode:raise RuntimeError('offset kernel/analytic/rehashed-checkpoint tests failed; see kernel.log')
    zero=zero_equivalence(output/'zero',binary)
    positive=study(output/'positive',binary,.05)
    negative=study(output/'negative',binary,-.05)
    refusal=probe(output/'retained-failures/scope-refusals',binary)
    result=dict(schema='soil-fem-plate6-curved-offset-acceptance/1',source_head='cecbb39edf7359f02ca3968998532ebd8612d544',
        passed=False,requested_all_items_passed=False,
        curved_tri6=dict(accepted=False,passed=False,reason='No validated covariant curved-shell kernel, nodal director basis or curved MITC extension; classical shell refinement benchmarks not run. Curved geometry remains explicitly refused.',benchmarks=refusal['classical_benchmarks']),
        thickness_offset=dict(accepted=True,passed=True,scope='Explicit mindlin_tri6_offset; direct straight/flat Tri6 input, finite |e|<=t/2, common signed normal offset at shared plate nodes; dry elastic MC soil + linear plate/beam; install/load/unload/remove only.',
            signed_distances_m=[-.05,.05],kernel_returncode=kernel.returncode,zero_equivalence=zero,
            studies={key:{k:v for k,v in value.items() if k!='reference'} for key,value in [('positive',positive),('negative',negative)]}),
        p6_complete=False,windows_ci='unconfirmed')
    write(output/'acceptance.json',result)
    metrics=dict(studies={k:json.loads((output/k/'acceptance.json').read_text()) for k in ('positive','negative')},zero=zero)
    write(output/'measurements.json',metrics)
    files=[]
    for p in sorted(output.rglob('*')):
        if p.is_file() and p.suffix!='.exe':files.append(dict(path=str(p.relative_to(output)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
    write(output/'evidence-manifest.json',dict(files=files,solver_sha256=hashlib.sha256(binary.read_bytes()).hexdigest()))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--binary',type=Path)
    args=p.parse_args();print(json.dumps(record(args.output,args.binary),indent=2))
