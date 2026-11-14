import sys,subprocess,time,json
from pathlib import Path
sys.path.insert(0,'interim/2026-12-19-c04-transient-water-cause-diagnosis')
from c04cause_driver import write
out=Path('work/c04-global-final');records=[]
for label in ('native-m128-w32','adaptive-m128-w32'):
 cmd=[sys.executable,'interim/2026-12-19-c04-transient-water-cause-diagnosis/c04cause_driver.py','--output',str(out),'--label',label,'--mechanical-steps','128','--water-steps','32']
 start=time.monotonic();peak=0
 with (out/(label+'.log')).open('x',encoding='utf-8',newline='\n') as log:
  child=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
  while child.poll() is None:
   try:
    for line in Path(f'/proc/{child.pid}/status').read_text().splitlines():
     if line.startswith(('VmHWM:','VmRSS:')):peak=max(peak,int(line.split()[1]))
   except FileNotFoundError:pass
   if time.monotonic()-start>900:child.kill();child.wait();break
   time.sleep(.1)
 records.append(dict(label=label,command=cmd,pid=child.pid,exit_code=child.returncode,elapsed_seconds=time.monotonic()-start,peak_sampled_rss_kib=peak))
 print(json.dumps(records[-1]),flush=True)
 if child.returncode:break
write(out/'m128-processes.json',dict(records=records,complete=len(records)==2 and all(r['exit_code']==0 for r in records),passed=False))
