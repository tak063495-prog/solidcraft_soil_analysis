"""Two curved attachment faces on the retained 12 Tet10, with soil prestrain."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import fem_reference as ref


def patch():
    mesh=json.loads((ROOT/'examples/plate6-offset-mesh.json').read_bytes())
    coords=np.array(mesh['nodes'])
    coords[:,2]+=.09*coords[:,0]**2-.05*coords[:,1]**2
    mesh['nodes']=coords.tolist()
    topology=ref.topology(mesh)
    owner={}
    faces=[]
    shared=[]
    for idx,tri in enumerate(mesh['face_sets']['top']):
        uid=f'f{idx}'
        cell=topology[tuple(sorted(tri))][0]
        ids=ref.face_nodes(mesh,cell,tri).tolist()
        xyz=coords[ids]
        if np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0])[2]<0:
            ids=[ids[k] for k in [0,2,1,5,4,3]]
        owner[uid]=dict(cell=int(cell),nodes=ids)
        for i in ids:
            if i not in shared:shared.append(i)
        faces.append(dict(id=uid,nodes=[shared.index(i) for i in ids],group='cap',offset_m=-.05))
    assert len(faces)==2
    base=6*len(shared)
    rows=[]
    for i in range(len(coords)):
        if i in shared: first=6*shared.index(i)
        else: first=base;base+=3
        rows.append(np.arange(first,first+3))
    rows=np.array(rows)
    mats={str(i):dict(young_kpa=i*12000.,poisson=.25) for i in mesh['material_ids']}
    q=np.zeros(base)
    q[rows[:,0]]=1e-4*coords[:,0]
    q[rows[:,2]]=-2e-4*coords[:,2]
    load=lambda t:dict(event='load',free_dofs=[],prescribed_configuration=q.tolist(),nodal_load=[0.]*base,plate_traction_kpa=[t,0.,0.])
    install=dict(event='install',holder='fixed soil',contacts={})
    remove=dict(event='remove',holder='fixed soil')
    steps=[load(0.),install,load(.5)]
    steps[-1]['free_dofs']=[int(rows[shared[3],0])]
    steps.extend([remove,install,load(-.25)])
    case=dict(nodes=coords[shared].tolist(),section=dict(id='liner',young_kpa=200000.,poisson=.25,thickness_m=.1,unit_weight_kn_m3=2.,drilling_factor=1e-4),
        faces=faces,directors=[],soil=dict(mesh=mesh,owners=owner,materials=mats,
        exposed_side={k:1 for k in owner},tie_at_physical_surface=True),steps=steps)
    return case,shared,rows


def verify(cli,out):
    out.mkdir(parents=True,exist_ok=False)
    case,shared,rows=patch()
    source=out/'case.json'
    source.write_text(json.dumps(case)+'\n',encoding='utf-8')
    cmd=[cli,'curved-shell',source]
    full=out/'full.json'
    subprocess.run([*cmd,full],check=True)
    frames=json.loads(full.read_bytes())['frames']
    mesh=case['soil']['mesh']
    max_gp=0.
    max_tie=0.
    free_residual=0.
    for frame in frames:
        disp=np.array(frame['soil_translations_m'])
        q=np.array(frame['configuration'])
        max_tie=max(max_tie,np.max(abs(disp[shared]-q[:6*len(shared)].reshape(-1,6)[:,:3])))
        for i,mat in enumerate(mesh['material_ids']):
            gp=ref.volume(mesh,i)
            c=ref.elastic_tensor(12000.*mat,.25)
            stress=np.einsum('ij,qjk,k->qi',c,gp['b'],disp[gp['nodes']].ravel())
            max_gp=max(max_gp,np.max(abs(stress-np.array(frame['soil_stress_kpa'][i]))))
    free_residual=abs(frames[3]['reactions'][case['steps'][2]['free_dofs'][0]])
    identical=[]
    for i in [2,3,4,5]:
        cp=out/f'cut-{i}.json'
        end=out/f'resume-{i}.json'
        subprocess.run([*cmd,cp,'--until',str(i)],check=True)
        subprocess.run([*cmd,end,'--resume',cp],check=True)
        identical.append(end.read_bytes()==full.read_bytes())
    result=dict(tet10=len(mesh['tetrahedra']),tri6=len(case['faces']),tie_max_m=float(max_tie),
        soil_gp_max_kpa=float(max_gp),free_residual_kn=free_residual,byte_equal=identical,
        frozen_equal=frames[3]['active']==frames[5]['archived'][0],generation=frames[5]['active']['generation'],
        initial_energy=sum(e['physical_energy_kn_m'] for e in frames[2]['active']['elements']),passed=False)
    (out/'measurements.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('cli',type=Path)
    p.add_argument('out',type=Path)
    args=p.parse_args()
    print(verify(args.cli.resolve(),args.out))
