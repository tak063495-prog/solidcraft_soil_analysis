"""Prepare, or explicitly recompute, the bounded274-cell independent benchmark."""
from pathlib import Path
import argparse,datetime,hashlib,importlib.metadata,json,os,shutil,subprocess,sys,time,zipfile
p=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--solver',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True,help='A new directory; existing evidence is never overwritten')
parser.add_argument('--run',action='store_true',help='Run both independent RK4 references and both Rust candidates; may take hours')
parser.add_argument('--increments',type=int,nargs=2,choices=[16,32,64],default=[32,64])
args=parser.parse_args();assert args.increments in [[16,32],[32,64]]
solver=args.solver.resolve();out=args.output.resolve()
read=lambda q:json.loads(q.read_text(encoding='utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
receipt=read(p/'preparation.json');manifest=read(p/'archive-manifest.json');archive=p/'immutable-inputs-and-reference.zip'
assert receipt['passed'] and sha(archive)==receipt['archive_sha256'] and sha(p/'archive-manifest.json')==receipt['manifest_sha256']
assert solver.is_file() and not out.exists() and out.parent.is_dir()
out.mkdir()
with zipfile.ZipFile(archive) as source:
 assert len(source.namelist())==len(set(source.namelist())) and set(source.namelist())==set(manifest)
 for n,expected in manifest.items():
  target=out/n;assert target.resolve().is_relative_to(out)
  data=source.read(n);assert len(data)==expected['bytes'] and hashlib.sha256(data).hexdigest()==expected['sha256']
  target.parent.mkdir(parents=True,exist_ok=True)
  with target.open('xb') as stream:stream.write(data)
for n in ['requirements.lock.txt','requirements-reference.lock.txt']:
 for pin in (out/'reference-code'/n).read_text(encoding='utf-8').splitlines():
  if '==' in pin and not pin.lstrip().startswith('#'):
   package,version=pin.strip().split('==');assert importlib.metadata.version(package)==version,(package,version)
tools=out/'reference-code/tools';inputs=out/'inputs'
sys.path.insert(0,str(tools))
from verify_mcc_structural_study import physical_case
base=read(inputs/'stages.json');bound={str(out/n):v['sha256'] for n,v in manifest.items()}
bound[str(solver)]=sha(solver);bound[str(Path(__file__).resolve())]=sha(Path(__file__))
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1')
records=[]
def call(label,command):
 command=list(map(str,command));started=time.monotonic()
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT)
  with (out/(label+'-child.json')).open('x',encoding='utf-8') as stream:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),stream)
  code=child.wait()
 records.append(dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-started));print(json.dumps(records[-1]),flush=True)
 assert code==0,label
for n in args.increments:
 folder=out/str(n);(folder/'ground').mkdir(parents=True)
 for name in ['mesh.json','ground.json']:shutil.copy2(inputs/'ground'/name,folder/'ground'/name)
 case=json.loads(json.dumps(base));case['settings']['initial_increment']=1/n
 assert physical_case(json.loads(json.dumps(case)))==physical_case(json.loads(json.dumps(base)))
 with (folder/'stages.json').open('x',encoding='utf-8') as stream:json.dump(case,stream,indent=2)
 for name in ['ground/mesh.json','ground/ground.json','stages.json']:bound[str(folder/name)]=sha(folder/name)
 call('prepare-native'+str(n),[solver,'prepare-stages-job',folder/'ground/mesh.json',folder/'ground/ground.json',folder/'stages.json',folder/'solver-job'])
if args.run:
 references=out/'references';references.mkdir()
 for n in [100,200]:call('reference'+str(n),[sys.executable,tools/'run_mcc_structural_reference.py',inputs/'ground/mesh.json',inputs/'ground/ground.json',inputs/'stages.json',references/('reference'+str(n)+'.json'),'--identity',inputs/'mesh-identity.json','--steps',n,'--integration-method','rk4','--checkpoint-dir',out/('reference'+str(n)+'-checkpoints')])
 for n in args.increments:
  folder=out/str(n);call('native'+str(n),[solver,'run-stages-job',folder/'solver-job']);shutil.copy2(folder/'solver-job/result.json',folder/'stage-results.json')
  call('export'+str(n),[sys.executable,tools/'export_stages.py',folder/'ground/mesh.json',folder/'stage-results.json',folder/'results'])
  call('audit'+str(n),[sys.executable,tools/'verify_stage_package.py',folder/'ground/mesh.json',folder/'ground/ground.json',folder/'stage-results.json',folder/'results','--solver',solver])
 call('independent-comparison',[sys.executable,tools/'verify_mcc_structural_study.py',out/str(args.increments[0]),out/str(args.increments[1]),out/'independent-verification.json','--references',references,'--reference-inputs',inputs])
assert all(sha(Path(n))==h for n,h in bound.items())
with (out/'reproduction.json').open('x',encoding='utf-8') as stream:json.dump(dict(prepared=True,increment_denominators=args.increments,executed=args.run,independent_comparison_passed=args.run and (out/'independent-verification.json').exists(),original_acceptance_thresholds_preserved=True,files_sha256=bound,records=records,archive_sha256=sha(archive),spatial_convergence_accepted=False,new_pc_acceptance=False,p6_complete=False),stream,indent=2)
print('PASS '+('full bounded independent benchmark' if args.run else 'preparation only; use --run with a different new output folder for full recomputation'),flush=True)
