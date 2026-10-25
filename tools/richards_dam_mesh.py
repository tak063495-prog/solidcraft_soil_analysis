"""Structured 3D unit-cube benchmark, with face-aligned reservoir edges."""
import argparse
from collections import defaultdict
import json
from pathlib import Path


def make(n,finite_width=False):
    root=Path(__file__).resolve().parents[1]
    m=json.loads((root/'examples/seepage-column-tet4-mesh.json').read_text(encoding='utf-8'));assert 4<=n<=24 and n%4==0
    points=[[x/n,y/n,z/n] for z in range(n+1) for y in range(n+1) for x in range(n+1)];idx=lambda x,y,z:x+(n+1)*(y+(n+1)*z);tets=[];ids=[]
    for z in range(n):
        for y in range(n):
            for x in range(n):
                c=[idx(x,y,z),idx(x+1,y,z),idx(x+1,y+1,z),idx(x,y+1,z),idx(x,y,z+1),idx(x+1,y,z+1),idx(x+1,y+1,z+1),idx(x,y+1,z+1)]
                for f in [(0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)]:tets.append([c[i] for i in f]);ids.append(1 if z<n/2 else 2)
    m.update(nodes=points,tetrahedra=tets,material_ids=ids,region_ids=ids,schema='soil-fem-mesh/2',excavation_sets={});m.pop('tetrahedra10',None);m.pop('quadrature_order',None)
    owners=defaultdict(list)
    for ei,t in enumerate(tets):
        for f in [(1,2,3),(0,3,2),(0,1,3),(0,2,1)]:owners[tuple(sorted(t[i] for i in f))].append(ei)
    reservoir=[];air=[]
    for face,o in owners.items():
        if len(o)!=1:continue
        p=[points[i] for i in face];xmin=all(v[0]==0. for v in p);xmax=all(v[0]==1. for v in p);top=all(v[2]==1. for v in p)
        inside_width=not finite_width or (min(v[1] for v in p)>=.25 and max(v[1] for v in p)<=.75)
        if xmin and max(v[2] for v in p)<=.5 and inside_width:reservoir.append(list(face))
        elif xmax or top or xmin:air.append(list(face))
    m['face_sets']={'reservoir':reservoir,'air':air};m['node_sets']={'all':list(range(len(points))),'bottom':[i for i,p in enumerate(points) if p[2]==0.],'top':[i for i,p in enumerate(points) if p[2]==1.]}
    return m


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--subdivisions',type=int,required=True);p.add_argument('--finite-width',action='store_true');a=p.parse_args()
    with a.output.open('x',encoding='utf-8') as f:json.dump(make(a.subdivisions,a.finite_width),f)
