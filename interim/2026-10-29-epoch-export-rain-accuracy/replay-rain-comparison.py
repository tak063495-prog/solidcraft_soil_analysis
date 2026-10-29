"""Replay the frozen comparison from audited caches; no raw nonlinear rerun."""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys
source=Path(__file__).resolve().parent
read=lambda q:json.loads(q.read_text('utf-8'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
inventory=read(source/'files-sha256.json')
assert all(sha(source/n)==h for n,h in inventory.items())
assert len(sys.argv)==2
output=Path(sys.argv[1]).resolve();assert not output.exists() and not output.is_relative_to(source) and not source.is_relative_to(output)
shutil.copytree(source,output)
refined='richards-local-rain-3d-refinement-20261029-v1'
command=[sys.executable,'-X','utf8',str(output/refined/'compare-caches.py')]
completed=subprocess.run(command,cwd=output,capture_output=True,text=True,timeout=600)
(output/'comparison.log').write_text(completed.stdout+completed.stderr,'utf-8')
assert completed.returncode==0
expected=read(output/refined/'expected-comparison.json');actual=read(output/refined/'comparison.json')
def canonical(report):
    report=dict(report);bound={}
    for name,digest in report['audited_input_sha256'].items():
        name=name.replace('\\','/')
        if 'richards-local-rain-3d' in name:
            name=name[name.index('richards-local-rain-3d'):]
        else:name=refined+'/'+name
        assert name not in bound;bound[name]=digest
    report['audited_input_sha256']=bound
    return report
assert canonical(actual)==canonical(expected)
assert inventory=={n:sha(source/n) for n in inventory}
receipt=dict(passed=True,all_comparison_fields_exact_except_absolute_path_prefixes=True,
    original_comparison_source_unchanged=True,checks=len(actual['checks']),
    all_predeclared_gates_passed=actual['passed'],raw_physical_audits_reexecuted=False,
    nonlinear_solves_reexecuted=False,source_packet_unchanged=True,p6_complete=False)
(output/'portable-verification.json').write_text(json.dumps(receipt,indent=2),'utf-8')
print(json.dumps(receipt))
