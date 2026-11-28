"""Repeat only the existing small 5 kPa NC case; failures never become acceptance."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from probe_nc_midpoint import inputs
from probe_low_preload_nc import sha


def run(binary, output):
    output.mkdir(parents=True,exist_ok=False)
    processes=[]; binary_sha=sha(binary)
    for n in (128,256):
        folder=output/f'input{n}';folder.mkdir()
        mesh,ground,case=inputs('backward_euler');case['stages'][2]['initial_increment']=1/n
        for name,value in zip(('mesh','ground','stages'),(mesh,ground,case)):
            (folder/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
        for mode in ('trace','direct'):
            command=[str(binary),mode,*[str(folder/(name+'.json')) for name in ('mesh','ground','stages')],str(output/f'{mode}{n}')]
            start=time.monotonic()
            with (output/f'{mode}{n}.log').open('x') as log:
                child=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=False)
            processes.append(dict(command=command,exit_code=child.returncode,elapsed_seconds=time.monotonic()-start))
            (output/'processes.json').write_text(json.dumps(processes,indent=2)+'\n',encoding='utf-8')
            if child.returncode:raise ValueError('diagnostic process failed; no acceptance')
    for command in ([sys.executable,str(HERE/'extract_oblique.py'),'--output',str(output/'prior')],
                    [sys.executable,str(HERE/'verify_restarts.py'),'--data',str(output),'--binary',str(binary),'--output',str(output/'restarts')],
                    [sys.executable,str(HERE/'measure.py'),'--data',str(output),'--output',str(output/'measurements.json')]):
        subprocess.run(command,check=True)
    if sha(binary)!=binary_sha:raise ValueError('binary changed during study')
    print(json.dumps(dict(passed=False,solver_binary_sha256=binary_sha,
                         measurements=str(output/'measurements.json'),completion_restart_verified=False),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.binary.resolve(),a.output.resolve())
