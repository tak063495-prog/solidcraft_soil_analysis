from pathlib import Path
import argparse,subprocess,time,json,sys,os
root=Path.cwd();study=root/'interim/2026-11-08-nc-compensated-linear-replay'
p=argparse.ArgumentParser();p.add_argument('mode',choices=['rust','python']);p.add_argument('--output',type=Path,required=True);p.add_argument('--probe',default='work/nc-linear/candidate/target/release/examples/nc_stop_probe'+('.exe' if os.name=='nt' else ''));p.add_argument('--prior',default='work/nc-linear/evidence/prior');p.add_argument('--traces',default='work/nc-linear/evidence/diagnostic');a=p.parse_args();mode=a.mode;out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir()
checks={
 'rust':[
  ('root-fmt',['cargo','fmt','--all','--check']),
  ('rust-linalg',['cargo','test','--locked','--release','--lib','linalg::tests']),
  ('rust-free-contact',['cargo','test','--locked','--release','--test','free_contact']),
  ('rust-mcc-midpoint',['cargo','test','--locked','--release','--test','mcc_midpoint']),
 ],
 'python':[
  ('python-midpoint',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_nc_midpoint_resume.py','-v']),
  ('python-diagnosis',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_low_preload_nc_diagnosis.py','-v']),
  ('python-trace',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_nc_newton_trace.py','-v']),
  ('python-linear-audit',[sys.executable,str(study/'test_linear_audit.py'),'-v']),
  ('python-diagnostic-guards',[sys.executable,'interim/2026-11-08-nc-midpoint-stop-diagnosis/test_diagnosis.py','--probe',a.probe,'--prior',a.prior,'--traces',a.traces,'-v']),
 ]}
records=[]
for name,args in checks[mode]:
 start=time.monotonic()
 with (out/'logs'/f'{name}.log').open('x') as f:p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,timeout=300)
 records.append(dict(name=name,command=args,exit_code=p.returncode,elapsed_seconds=time.monotonic()-start))
 (out/f'{mode}-check-processes.json').write_text(json.dumps(records,indent=2)+'\n')
 print(name,p.returncode,flush=True)
 if p.returncode:sys.exit(p.returncode)
