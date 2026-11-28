import json,pathlib,subprocess,time
root=pathlib.Path('/workspace/solidcraft_soil_analysis');out=root/'work/a04cause-validation';out.mkdir(parents=True,exist_ok=True)
for name,cmd in [('fmt',['cargo','fmt','--all','--check']),('clippy',['cargo','clippy','--locked','--all-targets','--','-D','warnings']),('build',['cargo','build','--locked','--release','--all-targets']),('rust-discovery',['cargo','test','--locked','--release','--all-targets','--no-run','--message-format=json']),('rust-main',['cargo','test','--locked','--release']),('rust-examples',['cargo','test','--locked','--release','--examples'])]:
 start=time.monotonic()
 with (out/(name+'.log')).open('x') as log:
  p=subprocess.Popen(cmd,cwd=root,stdout=log,stderr=subprocess.STDOUT);code=p.wait()
 row=dict(command=cmd,pid=p.pid,exit_code=code,elapsed_seconds=time.monotonic()-start)
 (out/(name+'-process.json')).write_text(json.dumps(row,indent=2)+'\n');print(name,row,flush=True)
