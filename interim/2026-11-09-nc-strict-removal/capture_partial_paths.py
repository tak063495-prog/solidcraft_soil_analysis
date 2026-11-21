"""Current-head partial step studies; no full-path or acceptance claim."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from probe_nc_midpoint import inputs
from probe_low_preload_nc import write, read, sha


def capture(out, solver):
    out.mkdir(parents=True,exist_ok=False)
    records=[]
    def run(args,label,timeout):
        start=time.monotonic(); args=list(map(str,args))
        with (out/(label+'.log')).open('x') as log:
            r=subprocess.run(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        records.append(dict(command=args,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start))
        (out/'processes.json').write_text(json.dumps(records,indent=2)+'\n');r.check_returncode()
    for count in (128,256):
        folder=out/('native-'+str(count)); folder.mkdir()
        values=inputs('midpoint');values[2]['stages'][2]['initial_increment']=1/count
        for n,v in zip(('mesh','ground','stages'),values):write(folder/(n+'.json'),v)
        for progress in (.875,.9375):
            # Existing verification rejection halves the first increment; the
            # actual checkpoint's stage/progress is checked, never assumed.
            accepted=14+int(progress*count*2)
            dest=folder/(str(progress)+'.json')
            run([solver,'checkpoint-stages',folder/'mesh.json',folder/'ground.json',folder/'stages.json',accepted,dest],
                'native-'+str(count)+'-'+str(progress),90)
            r=read(dest)['run'];assert r['last_committed']['stage_index']==3 and r['last_committed']['stage_progress']==progress
        write(folder/'input-sha256.json',{n+'.json':sha(folder/(n+'.json')) for n in ('mesh','ground','stages')})
    for count in (128,256):
        run([sys.executable,ROOT/'tools/low_preload_nc_reference.py','--steps',count,'--method','midpoint',
             '--budget-seconds',60,'--output',out/('reference-'+str(count))],'reference-'+str(count),75)
    write(out/'solver-sha256.json',{'sha256':sha(solver)})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--solver',type=Path,required=True)
    a=p.parse_args();capture(a.output.resolve(),a.solver.resolve())
