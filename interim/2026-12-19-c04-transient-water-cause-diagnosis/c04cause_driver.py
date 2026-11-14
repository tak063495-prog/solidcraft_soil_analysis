"""Opt-in small C-04 global controls; no production dispatch imports this module."""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time
import traceback

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OLD=HERE.parent/'2026-11-22-wet-joint-transient-water'
sys.path.insert(0,str(HERE))
import c04cause_reference as reference


def legacy():
    spec=importlib.util.spec_from_file_location('c04cause_legacy_study_helpers',OLD/'study.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def write(path,value):
    path=Path(path)
    with path.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,indent=2,ensure_ascii=False,allow_nan=False,default=convert)
        stream.write('\n')


def convert(value):
    import numpy as np
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def initialize(out,probe,solver,element='tet4'):
    s=legacy();out.mkdir(parents=True)
    mesh=s.read(OLD/'inputs/mesh.json') if element=='tet4' else s.read(ROOT/'examples/joint-water-column-tet10-mesh.json')
    write(out/'mesh.json',mesh);write(out/'ground.json',s.read(OLD/'inputs/ground.json'))
    for label,path in (('probe',probe),('solver',solver)):
        shutil.copy2(path,out/(label+path.suffix))
    import subprocess
    write(out/'identity.json',dict(starting_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
          element=element,binaries={p.name:s.sha(p) for p in out.iterdir() if p.name in ('probe','solver','probe.exe','solver.exe')},
          prescribed_quasistatic_water=True,hydraulic_joint_coupling=False,production_changed=False))


def solve(out,label,mechanical_steps,water_steps,parts,timeout,boost=False):
    s=legacy();folder=out/label
    if label.startswith('native'):
        s.native(out,folder,mechanical_steps,water_steps,timeout)
        return
    folder.mkdir()
    case=s.mechanical_case(mechanical_steps,water_steps)
    write(folder/'case.json',case)
    model=reference.Model(s.read(out/'mesh.json'),s.read(out/'ground.json'),case)
    start=time.monotonic()
    try:
        with (folder/'fields.jsonl').open('x',encoding='utf-8',newline='\n') as stream:
            def sink(row,fields):
                stream.write(json.dumps(dict(row=row,fields=fields),default=convert,allow_nan=False)+'\n')
                stream.flush()
            result=reference.run(model,mechanical_steps,case['settings']['material_maximum_substep'],
                                 case['joint_system']['maximum_increment_m'],sink,fixed_contact_parts=parts,boost_nontrivial=boost)
        write(folder/'result.json',result)
        write(folder/'process.json',dict(pid=os.getpid(),elapsed_seconds=time.monotonic()-start,
              native_results_used_as_input=False,global_equilibrium_solved=True,complete=True))
    except BaseException:
        write(folder/'failure.json',dict(complete=False,passed=False,traceback=traceback.format_exc(),
              elapsed_seconds=time.monotonic()-start))
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--initialize',action='store_true')
    p.add_argument('--probe',type=Path,default=ROOT/'target/release/examples/wet_joint_transient_probe')
    p.add_argument('--solver',type=Path,default=ROOT/'target/release/soil-fem-3d')
    p.add_argument('--element',choices=['tet4','tet10'],default='tet4')
    p.add_argument('--label');p.add_argument('--mechanical-steps',type=int,default=32)
    p.add_argument('--water-steps',type=int,default=32);p.add_argument('--fixed-contact-parts',type=int)
    p.add_argument('--timeout',type=int,default=900);p.add_argument('--boost-nontrivial',action='store_true')
    a=p.parse_args()
    if a.initialize:initialize(a.output.resolve(),a.probe.resolve(),a.solver.resolve(),a.element)
    else:
        if a.label is None:p.error('--label required for solve')
        solve(a.output.resolve(),a.label,a.mechanical_steps,a.water_steps,a.fixed_contact_parts,a.timeout,a.boost_nontrivial)

if __name__=='__main__':main()
