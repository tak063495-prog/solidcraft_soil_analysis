"""Existing inputs only. Shared IDs and duplicated attachment IDs are distinct."""
import copy
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
B02=HERE.parent/'2026-12-01-tri6-offset-composite-b02'
sys.path.insert(0,str(B02))
from cases import ROOT, inputs, read, write
import fem_reference as ref


def make(kind='kink-zero',level=0,steps=4):
    m,g,c,a=inputs(level,steps,0.,False)
    cfg=c['beam_system'];template=copy.deepcopy(cfg['plates'][0]);xyz=np.asarray(m['nodes'])
    topology=ref.topology(m);plates=[];nodes={}
    faces=[]
    for face,owners in topology.items():
        if len(owners)!=1 and kind!='height-step':continue
        top=np.all(xyz[list(face),2]==1.)
        side=np.all(xyz[list(face),0]==1.)
        middle=kind=='height-step' and np.all(xyz[list(face),2]==.5)
        if top or middle or (side and kind in ('kink-zero','kink-split')):
            ids=ref.face_nodes(m,owners[0],list(face)).tolist()
            normal=np.cross(xyz[ids[1]]-xyz[ids[0]],xyz[ids[2]]-xyz[ids[0]])
            axis=2 if top or middle else 0
            if normal[axis]<0:ids=[ids[j] for j in (0,2,1,5,4,3)]
            faces.append((ids,'top' if top else ('middle' if middle else 'side')))
    for i,(ids,plane) in enumerate(faces):
        label=plane if kind=='kink-split' else (('positive' if xyz[ids[:3]].mean(axis=0)[0]>xyz[ids[:3]].mean(axis=0)[1] else 'negative') if kind=='step-split' else 'shared')
        names=[f'{label}:n{j}' for j in ids]
        e=0. if kind=='kink-zero' else .05
        if kind=='step-split':e=.05 if xyz[ids[:3]].mean(axis=0)[0]>xyz[ids[:3]].mean(axis=0)[1] else -.05
        if kind=='height-step' and plane=='middle':e=-.025
        if kind=='normal-sign-pair' and i%2:
            names=[names[j] for j in (0,2,1,5,4,3)];e=-.05
        p=copy.deepcopy(template);p.update(id=f'p{i}',nodes=names,attachment_offset_m=e)
        plates.append(p)
        for j in ids:
            name=f'{label}:n{j}';sel=f'b03:{j}';m['node_sets'][sel]=[j]
            nodes[name]=dict(id=name,position_m=m['nodes'][j],soil_node_set=sel)
    cfg['nodes']=list(nodes.values());cfg['plates']=plates
    names=[p['id'] for p in plates]
    c['stages'][0]['members']=names;c['stages'][-1]['members']=names
    return m,g,c,a


def detached_fold(angle=90.,e1=.05,e2=.05):
    m,g,c,a=make('kink-zero')
    r=np.deg2rad(angle);points=[[0,0,1],[1,0,1],[0,1,1],[0,-np.cos(r),1+np.sin(r)]]
    ids=[[0,1,2],[1,0,3]]
    nodes={f'v{i}':dict(id=f'v{i}',position_m=p,match_soil_node=False) for i,p in enumerate(points)}
    plates=[];template=c['beam_system']['plates'][0]
    for k,tri in enumerate(ids):
        names=[f'v{i}' for i in tri]
        for i,j in ((0,1),(1,2),(2,0)):
            pair=tuple(sorted((tri[i],tri[j])));name=f'm{pair[0]}-{pair[1]}'
            nodes[name]=dict(id=name,position_m=((np.array(points[pair[0]])+points[pair[1]])/2).tolist(),match_soil_node=False)
            names.append(name)
        p=copy.deepcopy(template);p.update(id=f'fold{k}',nodes=names,attachment_offset_m=(e1,e2)[k]);plates.append(p)
    cfg=c['beam_system'];cfg['nodes']=list(nodes.values());cfg['plates']=plates
    c['stages'][0]['members']=[p['id'] for p in plates];c['stages'][-1]['members']=[p['id'] for p in plates]
    return m,g,c,a
