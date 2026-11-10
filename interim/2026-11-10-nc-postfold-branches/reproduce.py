"""Bounded small-case experiments, at most two independent processes at a time."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path
import subprocess
import sys
import time
from paths import HERE


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def run(args):
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    records=[]
    def invoke(specification):
        label,driver,extra=specification
        command=[sys.executable,str(HERE/driver),'--evidence',str(args.evidence),'--binary',str(args.binary),
                 '--output',str(out/label),*extra]
        started=time.monotonic()
        with (out/(label+'.log')).open('x',encoding='utf-8') as log:
            result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        record=dict(label=label,command=command,exit_code=result.returncode,elapsed_seconds=time.monotonic()-started)
        print(json.dumps(record),flush=True)
        return record
    def batch(specifications):
        with ThreadPoolExecutor(max_workers=2) as pool:
            for record in pool.map(invoke,specifications):
                records.append(record)
                (out/'processes.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
        if any(record['exit_code'] for record in records):
            raise ValueError('diagnostic subprocess failed; preserve its log')
    batch([
        ('native-long','trace.py',['--map','native','--step-size','.0002','--count','2048','--budget-seconds','600']),
        ('independent-long','trace.py',['--map','independent','--step-size','.0002','--count','768','--budget-seconds','600']),
        ('native-other-anchor','trace.py',['--map','native','--anchor','.97','--step-size','.0002','--count','512','--budget-seconds','240']),
    ])
    probes=[]
    for kind in ['native','independent']:
        folder=out/(kind+'-long')
        with gzip.open(folder/'samples.jsonl.gz','rt',encoding='utf-8') as stream:
            samples=[json.loads(line) for line in stream]
        if not samples:continue
        peak=max(samples,key=lambda p:p['coordinate'])
        for amplitude in ([-1.,-.25,.25,1.] if kind=='native' else [-.25,.25]):
            probes.append((f'{kind}-probe-{amplitude}','trace.py',['--map',kind,'--anchor-samples',str(folder/'samples.jsonl.gz'),
                '--anchor-step',str(peak['step']),'--probe',str(amplitude),'--step-size','.0001',
                '--count','64','--budget-seconds','120']))
    batch(probes)
    batch([(f'native-{order}','reorder.py',['--map','native','--order',order,'--budget-seconds','180'])
           for order in ['x-first','x-last','y-first','y-last','z-first','z-last']])
    # A Cartesian ordering has its own freshly integrated history, not a copied
    # state from the old complement-first path. Try arcs from the two highest s.
    candidates=sorted([(read(out/f'native-{order}'/'report.json')['final_coordinate'],order)
                       for order in ['x-first','x-last','y-first','y-last','z-first','z-last']],reverse=True)
    batch([(f'native-{order}-arc','trace.py',['--map','native','--order',order,
            '--anchor-samples',str(out/f'native-{order}'/'samples.jsonl.gz'),'--anchor','1.',
            '--step-size','.0002','--count','512','--budget-seconds','240']) for _,order in candidates[:2]])
    # Real middle checkpoint, same source/settings, separate child process.
    if (out/'native-long/middle.json').exists():
        batch([('native-long-resumed','trace.py',['--map','native','--step-size','.0002','--count','2048',
                '--budget-seconds','600','--resume',str(out/'native-long/middle.json')])])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['evidence','binary','output']:p.add_argument('--'+name,type=Path,required=True)
    run(p.parse_args())
