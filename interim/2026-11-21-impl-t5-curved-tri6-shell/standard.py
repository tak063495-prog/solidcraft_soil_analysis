"""Run the fixed B-01 structural reference cards through the production CLI."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

HERE = Path(__file__).resolve().parent
CARDS = HERE.parent / '2026-12-03-tri6-curved-b01/draft-reference/benchmarks.json'


def connectivity(n):
    def node(i, j):
        return (2*n+1)*j+i
    faces = []
    for j in range(0, 2*n, 2):
        for i in range(0, 2*n, 2):
            for offsets in [[(0,0),(2,0),(2,2),(1,0),(2,1),(1,1)],
                            [(0,0),(2,2),(0,2),(1,1),(1,2),(0,1)]]:
                ids = [node(i+a,j+b) for a,b in offsets]
                faces.append(dict(id=f'f{len(faces):04}', group='shell', nodes=ids, offset_m=0.))
    return faces


def build(card, n, distorted):
    name = card['id']
    mat, geom = card['material'], card['geometry']
    nodes, directors, fixed = [], [], set()
    for j in range(2*n+1):
        for i in range(2*n+1):
            k=len(nodes)
            u,v=i/(2*n),j/(2*n)
            if distorted and 0<i<2*n and 0<j<2*n:
                # Interior midsides belong to the same parametric edge on both faces.
                u+=0.10/n*np.sin(i*np.pi/2)*np.sin(j*np.pi/2)
            if name=='hemisphere_hole':
                phi=np.pi/2*u
                polar=np.pi/2-np.deg2rad(90-geom['hole_polar_angle_deg'])*v
                d=np.array([np.sin(polar)*np.cos(phi),np.sin(polar)*np.sin(phi),np.cos(polar)])
                x=geom['radius']*d
                if i==0: fixed.update(6*k+a for a in [1,3,5])
                if i==2*n: fixed.update(6*k+a for a in [0,4,5])
                if i==j==0: fixed.add(6*k+2)
            else:
                angle=np.deg2rad(geom.get('half_angle_deg',90))*u
                d=np.array([np.sin(angle),0.,np.cos(angle)])
                x=geom['radius']*d+[0,geom['length']*v/2,0]
                if i==0: fixed.update(6*k+a for a in [0,4,5])
                if i==2*n and name=='pinched_cylinder': fixed.update(6*k+a for a in [2,3,4])
                if j==0: fixed.update(6*k+a for a in [1,3,5])
                if j==2*n: fixed.update(6*k+a for a in [0,2,4])
            nodes.append(x.tolist())
            directors.append(dict(group='shell',node=k,direction=d.tolist()))
    faces=connectivity(n)
    num=6*len(nodes)
    load=np.zeros(num)
    traction=[0.,0.,0.]
    if name=='scordelis_lo': traction[2]=-90.
    if name=='pinched_cylinder': load[2]=-0.25
    if name=='hemisphere_hole':
        load[0]=1.
        load[6*(2*n)+1]=-1.
    free=sorted(set(range(num))-fixed)
    steps=[dict(event='install',contacts={},holder='diaphragm'),
           dict(event='load',prescribed_configuration=[0.]*num,free_dofs=free,
                nodal_load=load.tolist(),plate_traction_kpa=traction)]
    section=dict(id='benchmark',young_kpa=mat['E'],poisson=mat['nu'],
                 thickness_m=geom['thickness'],unit_weight_kn_m3=0.,drilling_factor=1e-4)
    return dict(nodes=nodes,faces=faces,directors=directors,section=section,steps=steps)


def execute(cli, out, case):
    source=out.with_suffix('.input.json')
    source.write_text(json.dumps(case,allow_nan=False)+'\n',encoding='utf-8')
    result=subprocess.run([cli,'curved-shell',source,out],capture_output=True,text=True)
    if result.returncode:
        return None,result.stderr.strip()
    return json.loads(out.read_bytes())['frames'][-1],None


def series(cli, folder):
    folder.mkdir(parents=True,exist_ok=False)
    rows=[]
    for card in json.loads(CARDS.read_bytes())['benchmarks']:
        for distorted in [False,True]:
            previous=None
            for n in [2,4,8]:
                case=build(card,n,distorted)
                xyz=np.array(case['nodes'])
                h=max(np.linalg.norm(xyz[a]-xyz[b]) for f in case['faces'] for a in f['nodes'] for b in f['nodes'])
                out=folder/f"{card['id']}-{distorted}-{n}.json"
                frame,error=execute(cli,out,case)
                row=dict(case=card['id'],n=n,distorted=distorted,h=float(h),diagnostic=error,
                         order=None,passed=False,reference=card['reference'])
                if frame:
                    probe=card['probe']
                    ids=np.flatnonzero(np.linalg.norm(xyz-np.array(probe['xyz']),axis=1)<probe['coordinate_tolerance'])
                    assert len(ids)==1
                    coordinate={'ux':0,'uy':1,'uz':2}[probe['component']]
                    u=frame['configuration'][6*ids[0]+coordinate]
                    ref=card['reference']['signed_value']
                    row.update(displacement=u,relative_difference=(u-ref)/abs(ref),
                               difference_from_previous=None if previous is None else u-previous,
                               physical_energy=sum(e['physical_energy_kn_m'] for e in frame['active']['elements']),
                               drilling_energy=sum(e['drilling_energy_kn_m'] for e in frame['active']['elements']),
                               external_work=frame['cumulative_applied_work_kn_m'],
                               free_residual=float(np.linalg.norm(np.array(frame['reactions'])[case['steps'][1]['free_dofs']])))
                    previous=u
                rows.append(row)
                (folder/'measurements.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
                print(card['id'],n,distorted,row.get('displacement'),error,flush=True)
    return rows


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('cli',type=Path)
    parser.add_argument('folder',type=Path)
    args=parser.parse_args()
    series(args.cli.resolve(),args.folder)
