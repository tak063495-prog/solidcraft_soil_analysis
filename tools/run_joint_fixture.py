"""Capture the small 3D Tet10 interface specimen in fresh CLI processes."""
import argparse
import json
from pathlib import Path
import shutil
import sys
from verify_checkpoint import invoke


def run(solver,ground,stages,output,mesh=None):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False);(output/'ground').mkdir()
    mp=output/'ground/mesh.json';gp=output/'ground/ground.json';sp=output/'stages.json';rp=output/'stage-results.json'
    shutil.copyfile(mesh or root/'examples/joint-column-tet10-mesh.json',mp);shutil.copyfile(ground,gp);shutil.copyfile(stages,sp)
    invoke([str(solver),'solve-stages',str(mp),str(gp),str(sp),str(rp)],output/'solve.log')
    invoke([sys.executable,str(root/'tools/export_stages.py'),str(mp),str(rp),str(output/'results')],output/'export.log')
    invoke([sys.executable,str(root/'tools/verify_stage_package.py'),str(mp),str(gp),str(rp),str(output/'results'),'--solver',str(solver)],output/'audit.log')
    native=json.loads(rp.read_text());print('JOINT fixture accepted:',output,'frames:',len(native['frames']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--ground',type=Path,required=True);p.add_argument('--stages',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--mesh',type=Path)
    a=p.parse_args();run(a.solver.resolve(),a.ground.resolve(),a.stages.resolve(),a.output.resolve(),a.mesh.resolve() if a.mesh else None)
