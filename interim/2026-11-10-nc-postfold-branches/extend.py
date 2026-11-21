"""Additional post-fold length and transverse probes on the computed z-last path."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent


def run(args):
    folder=args.study/'native-z-last-arc'
    with gzip.open(folder/'samples.jsonl.gz','rt',encoding='utf-8') as stream:
        samples=[json.loads(line) for line in stream]
    peak=max(samples,key=lambda p:p['coordinate'])
    specs=[('native-z-last-long',samples[-1]['step'],0.,2048,.0002,600)]
    specs += [(f'native-z-last-probe-{amplitude}',peak['step'],amplitude,128,.0001,180) for amplitude in [-.25,.25]]
    def invoke(spec):
        label,anchor,probe,count,step,budget=spec
        command=[sys.executable,str(HERE/'trace.py'),'--evidence',str(args.evidence),'--binary',str(args.binary),
                 '--output',str(args.study/label),'--map','native','--order','z-last',
                 '--anchor-samples',str(folder/'samples.jsonl.gz'),'--anchor-step',str(anchor),
                 '--probe',str(probe),'--count',str(count),'--step-size',str(step),'--budget-seconds',str(budget)]
        started=time.monotonic()
        with (args.study/(label+'.log')).open('x',encoding='utf-8') as log:
            result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        record=dict(label=label,command=command,exit_code=result.returncode,elapsed_seconds=time.monotonic()-started)
        print(json.dumps(record),flush=True)
        return record
    with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(invoke,specs))
    (args.study/'extension-processes.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
    if any(record['exit_code'] for record in records):raise ValueError('extension child failed; preserve log')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['evidence','binary','study']:p.add_argument('--'+name,type=Path,required=True)
    run(p.parse_args())
