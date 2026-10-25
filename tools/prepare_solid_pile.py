"""Capture and mesh the finite, perfectly bonded solid single-pile example.

This is a gravity-equilibrated elastic control, not a capacity calculation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tools'))
import fem_reference as ref

def write(path,value):
    with path.open('x',encoding='utf-8') as stream:json.dump(value,stream,indent=2,allow_nan=False)
def call(*args):
    subprocess.run(list(map(str,args)),cwd=root,check=True)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

parser=argparse.ArgumentParser()
parser.add_argument('output',type=Path)
parser.add_argument('--cli',type=Path,required=True)
parser.add_argument('--coarse-mesh-size-m',type=float,default=3.)
parser.add_argument('--fine-mesh-size-m',type=float,default=2.)
args=parser.parse_args()
import math
if not all(math.isfinite(v) for v in (args.fine_mesh_size_m,args.coarse_mesh_size_m)) or not 0 < args.fine_mesh_size_m < args.coarse_mesh_size_m:
    raise ValueError('require finite 0 < fine mesh size < coarse mesh size')
folder=args.output.resolve();folder.mkdir(parents=True,exist_ok=False)
for name,target in (('solid-pile.json','design.json'),('ground-solid-pile.json','ground.json'),('stages-solid-pile.json','stages.json')):
    shutil.copy2(root/'examples'/name,folder/target)
cli = args.cli.resolve()
cad = folder / 'cad'
call(sys.executable,'-X','utf8',root/'tools/cad_layers.py','--cli',cli,'--output-dir',cad,'--ground',folder/'ground.json','--design',folder/'design.json')
for size, label in ((args.coarse_mesh_size_m,'coarse'),(args.fine_mesh_size_m,'fine')):
    destination = folder / label
    destination.mkdir()
    mesh_path = destination / 'mesh.json'
    call(sys.executable,'-X','utf8',root/'tools/mesh_layers.py',cad/'bodies.json',mesh_path,'--mesh-size-m',size,'--order',2,'--quadrature-order',3)
    mesh = json.loads(mesh_path.read_text())
    topology = ref.topology(mesh)
    head, tip, surface = [], [], []
    for face, owners in topology.items():
        material = [mesh['material_ids'][e] for e in owners]
        xyz = [mesh['nodes'][i] for i in face]
        if len(owners) == 1 and material == [2] and all(abs(p[2]-8.) < 1e-10 for p in xyz):
            head.append(list(face))
        if sorted(material) == [1,2] and all(abs(p[2]-2.) < 1e-10 for p in xyz):
            tip.append(list(face))
        if len(owners) == 1 and material == [1] and all(abs(p[2]-8.) < 1e-10 for p in xyz):
            surface.append(list(face))
    assert head and tip and surface
    for name, faces in (('Pile:head',head),('Pile:tip',tip),('Soil:surface',surface)):
        mesh.setdefault('face_sets',{})[name] = faces
        ids = set()
        for face in faces:
            owner = topology[tuple(sorted(face))][0]
            ids.update(map(int,ref.face(mesh,owner,face)['nodes']))
        mesh['node_sets'][name] = sorted(ids)
    # Retain the exact raw mesher file and distinguish the added load/monitor sets.
    mesh_path.rename(destination / 'mesh-before-tags.json')
    write(mesh_path,mesh)
    write(destination / 'input-binding.json',dict(schema='soil-fem-solid-pile-input/1',
        cad_manifest_sha256=sha(cad/'bodies.json'),geometry_is_finite_3d=True,
        perfect_bond_shared_nodes=True,plane_strain_row_conversion=False,
        pile_dimensions_m=[1.,1.,6.],head_traction_kpa=[10.,0.,-100.],
        mesh_size_m=size,nodes=len(mesh['nodes']),tet10=len(mesh['tetrahedra']),
        files_sha256={str(p.relative_to(folder)):sha(p) for p in (mesh_path,folder/'ground.json',folder/'stages.json',folder/'design.json')},
        accepted=False,limitations='Gravity equilibrated elastic increment; no slip, capacity, construction effects, spacing conversion or spatial convergence acceptance'))
    print(label, len(mesh['nodes']), len(mesh['tetrahedra']), flush=True)
