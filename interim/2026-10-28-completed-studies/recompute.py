"""Recompute one fully specified dry MCC benchmark; no CAD mesher required."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,time,zipfile
p=Path(__file__).resolve().parent;read=lambda q:json.loads(q.read_text(encoding='utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
cases=read(p/'cases.json');parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--case',required=True,choices=sorted(cases));parser.add_argument('--solver',required=True,type=Path);parser.add_argument('--output',required=True,type=Path);parser.add_argument('--run',action='store_true',help='Run the nonlinear construction history; otherwise prepare only')
args=parser.parse_args();out=args.output.resolve();solver=args.solver.resolve();assert solver.is_file() and not out.exists() and out.parent.is_dir()
proof=read(p/'preparation.json');archive=p/'unsolved-inputs-and-native-models.zip';manifest=read(p/'archive-manifest.json')
assert proof['passed'] and sha(archive)==proof['archive_sha256'] and sha(p/'archive-manifest.json')==proof['manifest_sha256'] and sha(p/'cases.json')==proof['cases_sha256']
out.mkdir();inputs=out/'inputs';inputs.mkdir();names=['mesh.json','ground.json','stages.json']
with zipfile.ZipFile(archive) as z:
 assert len(z.namelist())==len(set(z.namelist())) and set(z.namelist())==set(manifest)
 for n in names:
  path='cases/'+args.case+'/'+n;data=z.read(path);assert len(data)==manifest[path]['bytes'] and hashlib.sha256(data).hexdigest()==manifest[path]['sha256']
  with (inputs/n).open('xb') as stream:stream.write(data)
binding={n:sha(inputs/n) for n in names};solver_sha=sha(solver);records=[]
def call(label,command):
 start=time.monotonic()
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:result=subprocess.run(list(map(str,command)),stdout=log,stderr=subprocess.STDOUT)
 records.append(dict(step=label,exit_code=result.returncode,elapsed_seconds=time.monotonic()-start));assert result.returncode==0,label
job=out/'solver-job';call('prepare',[solver,'prepare-stages-job',inputs/'mesh.json',inputs/'ground.json',inputs/'stages.json',job])
if args.run:
 call('solve',[job/'solver/soil-fem-3d.exe','run-stages-job',job]);shutil.copy2(job/'result.json',out/'stage-results.json')
 result=read(out/'stage-results.json');assert result['complete'] and len(result['frames'])==5
assert all(sha(inputs/n)==h for n,h in binding.items()) and sha(solver)==solver_sha
with (out/'recomputation.json').open('x',encoding='utf-8') as stream:json.dump(dict(prepared=True,executed=args.run,case=args.case,original_input_bytes_exact=True,files_sha256=binding,solver_sha256=solver_sha,reference_solver_sha256=cases[args.case]['solver_binary_sha256'],same_binary_as_reference=solver_sha==cases[args.case]['solver_binary_sha256'],records=records,requires_further_independent_field_audits_after_recomputation=True,spatial_convergence_accepted=False,new_pc_acceptance=False,p6_complete=False),stream,indent=2)
print('PASS '+('nonlinear solve complete; inspect/audit the new result' if args.run else 'preparation only; use --run with a different new output folder to recompute'))
