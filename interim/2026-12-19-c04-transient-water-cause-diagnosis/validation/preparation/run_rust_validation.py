import json,subprocess,time,sys
from pathlib import Path
out=Path('work/c04-software-validation');records=[]
jobs=[('fmt',['cargo','fmt','--all','--check']),('clippy',['cargo','clippy','--locked','--all-targets','--','-D','warnings']),('build',['cargo','build','--locked','--release']),('rust-list',['cargo','test','--locked','--release','--','--list']),('rust-all',['cargo','test','--locked','--release'])]
for name,cmd in jobs:
 start=time.monotonic()
 with (out/(name+'.log')).open('x',encoding='utf-8',newline='\n') as log:
  child=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT);rc=child.wait()
 records.append(dict(name=name,command=cmd,pid=child.pid,exit_code=rc,elapsed_seconds=time.monotonic()-start))
 (out/'rust-processes.json').write_text(json.dumps(records,indent=2)+'\n')
 print(name,rc,flush=True)
 if rc:sys.exit(rc)
