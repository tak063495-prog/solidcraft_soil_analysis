from pathlib import Path
import sys,json,zipfile,hashlib,subprocess,shutil
import argparse,os
parser=argparse.ArgumentParser();parser.add_argument('--baseline',type=Path,required=True);parser.add_argument('--current',type=Path);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
root=Path(__file__).resolve().parents[2];out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
current=args.current or root/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')
solvers={k:out/(k+'-solver'+('.exe' if os.name=='nt' else '')) for k in ('baseline','current')}
for k,source in [('baseline',args.baseline),('current',current)]:shutil.copy2(source,solvers[k])
def semantic(v):
 if isinstance(v,dict):return {k:semantic(x) for k,x in v.items() if k not in ('elapsed_seconds','solver_binary_sha256')}
 if isinstance(v,list):return [semantic(x) for x in v]
 return v
cases=[]
z=zipfile.ZipFile(root/'interim/2026-11-13-joint-ext-combos/legacy-evidence.zip')
for name in ('dry','wet','tip','mindlin-tri6'):cases.append(('legacy-'+name,{f:json.loads(z.read('baseline/'+name+'/'+f+'.json')) for f in ('mesh','ground','case')}))
z=zipfile.ZipFile(root/'interim/2026-11-13-joint-ext-combos/evidence.zip')
for name in ('tip-elastic-stress_free','tip-elastic-transfer_force','tip-history-stress_free','tip-history-transfer_force','wet-elastic','wet-hydrostatic'):cases.append(('scope-a-'+name,{f:json.loads(z.read(name+'/'+f+'.json')) for f in ('mesh','ground','case')}))
cases.append(('legacy-cst-dkt',dict(mesh=json.loads((root/'examples/hardening-column-tet10-mesh.json').read_text()),ground=json.loads((root/'examples/ground-plate-beam-column.json').read_text()),case=json.loads((root/'examples/stages-plate-beam-column.json').read_text()))))
records=[]
for name,data in cases:
 p=out/name;p.mkdir()
 for k,v in data.items():(p/(k+'.json')).write_text(json.dumps(v))
 runs={}
 for variant,binary in solvers.items():
  result=p/(variant+'.json');q=subprocess.run([str(binary),'solve-stages',*(str(p/(f+'.json')) for f in ('mesh','ground','case')),str(result)],capture_output=True,text=True)
  (p/(variant+'.log')).write_text(q.stdout+q.stderr);assert q.returncode==0,(name,variant,q.stderr)
  runs[variant]=json.loads(result.read_text());assert runs[variant]['complete'],name
 assert json.dumps(semantic(runs['baseline']),sort_keys=True,allow_nan=False)==json.dumps(semantic(runs['current']),sort_keys=True,allow_nan=False),name
 records.append(dict(case=name,passed=True,complete_runs_exact_except_elapsed_seconds_and_binary_sha=True,input_sha256={k:hashlib.sha256((p/(k+'.json')).read_bytes()).hexdigest() for k in data}))
 print(name,'exact',flush=True)
(out/'legacy-comparison.json').write_text(json.dumps(dict(passed=True,baseline_head='cecbb39edf7359f02ca3968998532ebd8612d544',solver_sha256={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in solvers.items()},cases=records),indent=2))
