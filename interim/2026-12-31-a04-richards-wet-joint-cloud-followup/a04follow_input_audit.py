"""Read-only physical-input equivalence under explicit E/F/Q permutations."""
import argparse
import copy
import hashlib
import itertools
from a04follow_common import *

def normalized(mesh,ground,case):
    m=copy.deepcopy(mesh);c=copy.deepcopy(case)
    cells=m['tetrahedra'];quadratic=m.get('tetrahedra10') or [None]*len(cells)
    physical=[dict(tet4=t,tet10=q,material=mat,region=reg) for t,q,mat,reg in zip(cells,quadratic,m['material_ids'],m['region_ids'])]
    for joint in c.get('joint_system',{}).get('joints',[]):
        for face in (joint.get('crack_tip') or {}).get('faces',[]):
            for key in ('owner_a','owner_b'):face[key]=physical[face[key]]
    for key in ('tetrahedra','tetrahedra10','material_ids','region_ids'):m.pop(key,None)
    m['physical_cells']=sorted(physical,key=lambda v:tuple(v['tet4']))
    m['face_sets']={k:sorted(v) for k,v in m['face_sets'].items()}
    return dict(mesh=m,ground=ground,case=c)

def audit(run):
    result={}
    for kind in ('soil','dry','hydro'):
        variants={}
        for e,f,q in itertools.product((0,1),repeat=3):
            name=f'e{e}f{f}q{q}';folder=run/'order'/kind/name
            paths=[folder/(k+'.json') for k in ('mesh','ground','case')]
            value=normalized(*map(read,paths));raw=json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()
            variants[name]=dict(physical_input_sha256=hashlib.sha256(raw).hexdigest(),raw_sha256={p.name:sha(p) for p in paths})
        hashes={v['physical_input_sha256'] for v in variants.values()}
        if len(hashes)!=1:raise ValueError(('physical input differs',kind,hashes))
        result[kind]=dict(all_physical_inputs_exact=True,explicit_tip_owner_aliases_remapped=True,variants=variants)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();write(a.output,audit(a.run))
