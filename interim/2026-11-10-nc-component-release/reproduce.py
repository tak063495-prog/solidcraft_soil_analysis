"""List or run the bounded small-case component-release diagnostics."""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

HERE=Path(__file__).resolve().parent


def commands(w):
    result=[]
    binary=w/'observer-target/release/examples'/('material_diagnostic.exe' if sys.platform=='win32' else 'material_diagnostic')
    for kind in ['native','independent']:
        for order in ['horizontal','vertical']:
            for step in [128,256]:
                label=f'{kind}-{order}{step}'
                cmd=[sys.executable,str(HERE/'study.py'),'--map',kind,'--order',order+'-first',
                     '--evidence',str(w/'origin'),'--binary',str(binary),'--reference',str(w/'reference400'),
                     '--output',str(w/label),'--step-size',str(1/step),'--budget-seconds',str(300 if step==128 else 600)]
                result.append((label,cmd))
        label=kind+'-horizontal-resumed'
        cmd=[sys.executable,str(HERE/'study.py'),'--map',kind,'--order','horizontal-first',
             '--evidence',str(w/'origin'),'--binary',str(binary),'--reference',str(w/'reference400'),
             '--output',str(w/label),'--resume',str(w/(kind+'-horizontal128')/'midpoint.json'),
             '--count','64','--budget-seconds','300']
        result.append((label,cmd))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--list-only',action='store_true')
    a=p.parse_args();w=a.output.resolve();cmds=commands(w)
    if a.list_only:print(json.dumps(cmds,indent=2))
    else:
        subprocess.run([sys.executable,str(HERE/'prepare.py'),'--output',str(w)],check=True)
        records=[]
        for label,cmd in cmds:
            start=time.monotonic()
            with (w/(label+'.log')).open('x') as log:
                r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
            records.append(dict(label=label,command=cmd,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start))
            (w/'numerical-processes.json').write_text(json.dumps(records,indent=2)+'\n')
            r.check_returncode()
        subprocess.run([sys.executable,str(HERE/'measure.py'),'--evidence',str(w),'--output',str(w/'measurements.json')],check=True)
