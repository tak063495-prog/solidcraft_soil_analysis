"""Existing accepted dry/static hydro controls; permutations only."""
import argparse
import sys
from contextlib import contextmanager
from a04follow_common import *
sys.path.insert(0,str(ROOT/'tools'))

def binding(kind):
    old=ROOT/'interim'/('2026-11-25-joint-birth-projection-a02' if kind=='dry' else '2026-11-26-wet-soil-birth-darcy-a03')
    sys.path.insert(0,str(old))
    s=load('a04follow_old_'+kind+'_study',old/'study.py')
    r=sys.modules['reference'];fields=sys.modules.get('fields')
    return s,r,fields

@contextmanager
def controls(reverse):
    import fem_reference as fem
    oldface=fem.face;oldjacobi=fem.jacobi
    def face(*a,**kw):
        q=oldface(*a,**kw)
        return {k:(v[::-1].copy() if k!='nodes' else v) for k,v in q.items()}
    def jacobi(n,alpha):
        x,w=oldjacobi(n,alpha)
        return (x[::-1].copy(),w[::-1].copy()) if n==4 and alpha in (0,1) else (x,w)
    if reverse:fem.face=face;fem.jacobi=jacobi
    try:yield
    finally:fem.face=oldface;fem.jacobi=oldjacobi

def layout(model):
    s=model.system
    ek=[sorted(c) for c in model.mesh.get('tetrahedra10',[]) or model.mesh['tetrahedra']]
    return dict(element_keys=ek,owner_elements=[[ek[p['owner_a']],ek[p['owner_b']]] for p in s.points],positions=[p['position_m'] for p in s.points],
                contact_keys=[(sorted(p['nodes_a']),sorted(p['nodes_b']),list(p['shape'])) for p in s.points],
                normal=s.normal,weight=s.weight,c=s.c)

def main():
    p=argparse.ArgumentParser();p.add_argument('--folder',type=Path,required=True);p.add_argument('--kind',choices=('dry','hydro'),required=True)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--engine',choices=('native','independent'),required=True);p.add_argument('--e',action='store_true');p.add_argument('--f',action='store_true');p.add_argument('--q',action='store_true')
    a=p.parse_args();s,r,fields=binding(a.kind);folder=a.folder
    if not folder.exists():
        folder.mkdir(parents=True)
        mesh,g,c=s.fixture('stress_free',4,True) if a.kind=='dry' else s.fixture(4,2.)
        if a.e:
            for k in ('tetrahedra','tetrahedra10','material_ids','region_ids'):
                if k in mesh:mesh[k].reverse()
            # Explicit owner indices are aliases of the same physical cells.
            for joint in c['joint_system']['joints']:
                for face in (joint.get('crack_tip') or {}).get('faces',[]):
                    for key in ('owner_a','owner_b'):face[key]=len(mesh['tetrahedra'])-1-face[key]
        if a.f:mesh['face_sets']={k:list(reversed(v)) for k,v in mesh['face_sets'].items()}
        for k,v in zip(('mesh','ground','case'),(mesh,g,c)):write(folder/(k+'.json'),v)
    with controls(a.q):
        m=r.Model(*[read(folder/(k+'.json')) for k in ('mesh','ground','case')])
        write(folder/('a04follow-'+a.engine+'-layout.json'),layout(m))
        if a.engine=='native':
            if a.kind=='dry':s.native(folder,a.probe,4)
            else:s.native(folder,a.probe,a.solver,4)
        else:
            if a.kind=='dry':
                with (folder/'independent-fields.jsonl').open('x',encoding='utf-8',newline='\n') as stream:
                    result=s.run(m,4,s.Sink(stream))
                write(folder/'independent-result.json',result)
            else:s.independent(folder,4)

if __name__=='__main__':main()
