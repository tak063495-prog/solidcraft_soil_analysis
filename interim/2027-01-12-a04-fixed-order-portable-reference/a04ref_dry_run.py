"""Re-solve the eight original dry permutations and capture open/.25 iterates."""
import argparse
import itertools
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
FOLLOW = HERE.parent/'2026-12-31-a04-richards-wet-joint-cloud-followup'


def run(out,probe):
    out.mkdir(parents=True,exist_ok=False)
    for e,f,q in itertools.product((0,1),repeat=3):
        folder = out/f'e{e}f{f}q{q}'
        flags = [flag for flag,enabled in zip(('--e','--f','--q'),(e,f,q)) if enabled]
        subprocess.run([sys.executable,FOLLOW/'a04follow_joint_order.py','--kind','dry',
            '--folder',folder,'--engine','independent','--probe',probe,'--solver',probe,*flags],check=True)
        env = dict(os.environ,A04REF_DRY_TRACE=str(folder/'iterates.jsonl'))
        env.pop('A04FOLLOW_Q',None)
        if q:
            env['A04FOLLOW_Q']='1'
        subprocess.run([probe,*[folder/(k+'.json') for k in ('mesh','ground','case')],
                        folder/'native','--whole'],env=env,check=True)
        print(folder.name,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--probe',required=True,type=Path)
    a=p.parse_args()
    run(a.output.resolve(),a.probe.resolve())
