"""Fresh small soil solves and input-only reference, with directional h_z."""
import argparse
import copy
from contextlib import contextmanager
import itertools
from functools import lru_cache
import sys
from a04follow_common import *
OLD=ROOT/'interim/2026-11-29-richards-wet-joint-a04'
sys.path.insert(0,str(ROOT/'tools'))
import a04follow_soil_reference as reference

def legacy():
    previous=sys.modules.get('soil_reference'); paths=sys.path[:]
    try:
        sys.modules['soil_reference']=reference
        sys.path.insert(0,str(OLD))
        return load('a04follow_old_soil_study',OLD/'study.py')
    finally:
        sys.path[:]=paths
        if previous is None:sys.modules.pop('soil_reference',None)
        else:sys.modules['soil_reference']=previous

def fixture(layers,n,e=False,f=False):
    s=legacy();mesh,ground,case=s.fixture(n,False)
    template=copy.deepcopy(mesh)
    mesh['nodes']=[[x,y,j/layers] for j in range(layers+1) for y in (0.,1.) for x in (0.,1.)]
    pattern=template['tetrahedra'][:6]
    mesh['tetrahedra']=[[int(v)+4*j for v in cell] for j in range(layers) for cell in pattern]
    for key in ('material_ids','region_ids'):
        mesh[key]=[template[key][0 if j<layers//2 else 6] for j in range(layers) for _ in range(6)]
    mesh['node_sets']={}
    xyz=np.asarray(mesh['nodes'])
    planes={'xmin':(0,0.),'xmax':(0,1.),'ymin':(1,0.),'ymax':(1,1.),'bottom':(2,0.),'top':(2,1.),'floor':(2,.5)}
    for k,(axis,value) in planes.items():mesh['node_sets'][k]=np.flatnonzero(xyz[:,axis]==value).tolist()
    mesh['node_sets']['all']=list(range(len(xyz)))
    keys=set(tuple(sorted(set(c)-{v})) for c in mesh['tetrahedra'] for v in c)
    mesh['face_sets']={k:[list(face) for face in sorted(keys) if all(xyz[v,axis]==value for v in face)] for k,(axis,value) in planes.items()}
    if e:
        for key in ('tetrahedra','material_ids','region_ids'):mesh[key].reverse()
    if f:
        mesh['face_sets']={k:list(reversed(v)) for k,v in mesh['face_sets'].items()}
    return mesh,ground,case

@contextmanager
def controls(order,reverse):
    original=reference.hyd.rule
    def rule(_):
        bary,w=original(order)
        return (bary[::-1].copy(),w[::-1].copy()) if reverse else (bary,w)
    original_local=reference.hyd.Model.local
    @lru_cache(maxsize=256)
    def cached(model,values,dtype):
        return original_local(model,np.asarray(values,dtype=dtype))
    def local(model,h):
        a=np.asarray(h)
        return cached(model,tuple(a),a.dtype.str)
    reference.hyd.rule=rule
    reference.hyd.Model.local=local
    try:yield
    finally:
        reference.hyd.rule=original
        reference.hyd.Model.local=original_local
        cached.cache_clear()

def layout(model):
    def h(m):
        return dict(elements=m.elements,positions=m.position,weights=m.weight,face_keys=m.keys,
                    owners=[[(int(m.elements[e]),int(j)) for e,j in owners] for owners in m.owners])
    return dict(element_keys=[sorted(c) for c in model.cells.tolist()],full=h(model.full),cut=h(model.cut))

def solve(folder,engine,order=3,reverse=False):
    mesh,ground,case=[read(folder/(k+'.json')) for k in ('mesh','ground','case')]
    with controls(order,reverse):
        model=reference.Model(mesh,ground,case)
        write(folder/(engine+'-layout.json'),layout(model))
        if engine=='independent':
            f,emit=sink(folder/'independent-fields.jsonl')
            try:r=model.solve(emit)
            finally:f.close()
            write(folder/'independent-result.json',r)
        else:
            s=legacy()
            # Saved production Fluid is the evidence; Model only reassembles skeleton observation.
            s.native_fields(folder)

def main():
    p=argparse.ArgumentParser();p.add_argument('--folder',type=Path,required=True);p.add_argument('--engine',choices=('native','independent'),required=True)
    p.add_argument('--order',type=int,choices=(3,5,7),default=3);p.add_argument('--reverse',action='store_true')
    a=p.parse_args();solve(a.folder,a.engine,a.order,a.reverse)

if __name__=='__main__':main()
