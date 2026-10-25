"""Capture a SolidCraft volume and attach an explicit flat-plate wall trace."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import plate_surfaces

root=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def prepare(output,cli,sizes):
    if len(sizes)!=2 or not all(math.isfinite(v) and v>0 for v in sizes) or sizes[0]<=sizes[1]:
        raise ValueError('require finite coarse > fine >0 mesh sizes')
    files=[('plate-wall.json','design.json'),('ground-plate-wall.json','ground.json'),
        ('stages-plate-wall.json','stages-template.json'),('plate-wall-surfaces.json','surface-spec.json')]
    cli=cli.resolve();output=output.resolve()
    if not cli.is_file():raise ValueError('SolidCraft CLI executable is missing')
    for source,_ in files:read(root/'examples'/source)
    output.mkdir(parents=True,exist_ok=False)
    for source,target in files:
        shutil.copyfile(root/'examples'/source,output/target)
    def call(*args):subprocess.run(list(map(str,args)),cwd=root,check=True)
    call(sys.executable,'-X','utf8',root/'tools/cad_layers.py','--cli',cli,'--output-dir',output/'cad',
        '--ground',output/'ground.json','--design',output/'design.json')
    for label,size in zip(('coarse','fine'),sizes,strict=True):
        run=output/label;run.mkdir();(run/'ground').mkdir()
        raw=run/'ground/mesh-before-plates.json'
        call(sys.executable,'-X','utf8',root/'tools/mesh_layers.py',output/'cad/bodies.json',raw,
            '--mesh-size-m',size,'--order',2,'--quadrature-order',3)
        mesh,case,binding=plate_surfaces.attach(read(raw),read(output/'stages-template.json'),read(output/'surface-spec.json'))
        write(run/'ground/mesh.json',mesh);write(run/'stages.json',case)
        shutil.copyfile(output/'ground.json',run/'ground/ground.json')
        binding.update(mesh_size_m=size,cli_sha256=sha(cli),cad_manifest_sha256=sha(output/'cad/bodies.json'),
            original_input_sha256={name:sha(output/name) for name in ('design.json','ground.json','stages-template.json','surface-spec.json')},
            source_sha256={str(p.resolve()):sha(p) for p in (Path(__file__),Path(plate_surfaces.__file__),Path(plate_surfaces.ref.__file__))},
            files_sha256={name:sha(run/name) for name in ('ground/mesh.json','ground/ground.json','ground/mesh-before-plates.json','stages.json')},
            nodes=len(mesh['nodes']),tet10=len(mesh['tetrahedra']),equilibrium_verified=False,spatial_convergence_accepted=False)
        write(run/'plate-surface-binding.json',binding)
        print(json.dumps(dict(mesh=label,size_m=size,nodes=len(mesh['nodes']),tet10=len(mesh['tetrahedra']),surfaces=binding['surfaces'])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);parser.add_argument('--cli',type=Path,required=True)
    parser.add_argument('--mesh-sizes-m',type=float,nargs=2,default=[4.,3.]);args=parser.parse_args()
    prepare(args.output,args.cli,args.mesh_sizes_m)
