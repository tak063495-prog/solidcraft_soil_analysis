"""Capture and audit an IDH staged fixture with immutable input copies."""
import argparse
from pathlib import Path
import shutil
import sys
from verify_checkpoint import invoke


def run(solver,mesh,ground,stages,output):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False);(output/'ground').mkdir()
    mp=output/'ground/mesh.json';gp=output/'ground/ground.json';sp=output/'stages.json';rp=output/'stage-results.json'
    for source,target in ((mesh,mp),(ground,gp),(stages,sp)):shutil.copyfile(source,target)
    invoke([str(solver),'solve-stages',str(mp),str(gp),str(sp),str(rp)],output/'solve.log')
    invoke([sys.executable,str(root/'tools/export_stages.py'),str(mp),str(rp),str(output/'results')],output/'export.log')
    invoke([sys.executable,str(root/'tools/verify_stage_package.py'),str(mp),str(gp),str(rp),str(output/'results'),'--solver',str(solver)],output/'audit.log')
    print('IDH staged fixture accepted:',output,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('solver','mesh','ground','stages','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();run(*[getattr(a,n).resolve() for n in ('solver','mesh','ground','stages','output')])
