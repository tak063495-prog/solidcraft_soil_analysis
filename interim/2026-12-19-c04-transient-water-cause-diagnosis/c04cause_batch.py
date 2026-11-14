"""Run the bounded diagnostic matrix in separate processes, retaining every failure."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
from c04cause_driver import write
HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--timeout',type=int,default=900);a=p.parse_args();out=a.output.resolve()
    jobs=[('native-w32',32,32,None),('adaptive-w32',32,32,None),('fixed1-w32',32,32,1),('fixed2-w32',32,32,2),
          ('fixed4-w32',32,32,4),('fixed8-w32',32,32,8),('native-m64-w32',64,32,None),
          ('adaptive-m64-w32',64,32,None),('fixed4-m64-w32',64,32,4),('boosted-w32',32,32,None)]
    records=[]
    for name,m,w,c in jobs:
        cmd=[sys.executable,str(HERE/'c04cause_driver.py'),'--output',str(out),'--label',name,
             '--mechanical-steps',str(m),'--water-steps',str(w),'--timeout',str(a.timeout)]
        if c is not None:cmd+=['--fixed-contact-parts',str(c)]
        if name.startswith('boosted'):cmd+=['--boost-nontrivial']
        start=time.monotonic();peak=0;timeout=False
        with (out/(name+'.log')).open('x',encoding='utf-8',newline='\n') as log:
            log.write(repr(cmd)+'\n');log.flush()
            child=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
            while child.poll() is None:
                status=Path(f'/proc/{child.pid}/status')
                try:
                    for line in status.read_text().splitlines():
                        if line.startswith(('VmHWM:','VmRSS:')):peak=max(peak,int(line.split()[1]))
                except (FileNotFoundError,PermissionError):pass
                if time.monotonic()-start>a.timeout:
                    timeout=True;child.kill();child.wait();break
                time.sleep(.1)
        record=dict(label=name,command=cmd,pid=child.pid,exit_code=child.returncode,timeout=timeout,
                    elapsed_seconds=time.monotonic()-start,peak_sampled_rss_kib=peak or None,
                    mechanical_steps=m,water_steps=w,fixed_contact_parts=c)
        write(out/(name+'-process.json'),record);records.append(record)
        print(json.dumps(record),flush=True)
        if child.returncode!=0:
            write(out/'matrix-failure.json',dict(complete=False,passed=False,records=records));return 1
    write(out/'matrix-processes.json',dict(complete=True,passed=False,global_controls=records))
    return 0

if __name__=='__main__':sys.exit(main())
