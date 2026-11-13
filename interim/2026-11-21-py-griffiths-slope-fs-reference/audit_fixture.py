"""Read-only verification of pinned input data; does not execute original solver."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def audit(path):
    c=json.loads((Path(__file__).parent/'case.json').read_text('utf-8'))
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=c['source']['sha256']: raise ValueError('fixture digest differs')
    f=next(f for f in json.loads(raw)['fixtures'] if f['name']==c['source']['fixture'])
    n=np.asarray(f['nodes'],float);xy=n[:,1:3]
    if len(n)!=661 or len(f['elements'])!=200 or f['materials'][0][:8]!=[1,100000,.3,1,20,20,10,0]:
        raise ValueError('material/mesh mismatch')
    # Q8 edges: each exterior edge must lie on one of the four verified lines.
    from collections import Counter
    edges=Counter(tuple(sorted(pair)) for e in f['elements'] for pair in zip(e[1:5],e[2:5]+e[1:2]))
    boundary=[]
    for (a,b),count in edges.items():
        if count==1:
            p=xy[[a-1,b-1]]
            if not any(np.allclose(v,0,atol=1e-10,rtol=0) for v in [p[:,0],p[:,1],p[:,1]-10,p[:,0]+2*p[:,1]-32]):
                raise ValueError('unverified exterior geometry')
            boundary.append([a,b])
    fixed_x=np.isclose(xy[:,0],0)|np.isclose(xy[:,1],0)
    fixed_y=np.isclose(xy[:,1],0)
    if not np.array_equal(n[:,3].astype(bool),fixed_x) or not np.array_equal(n[:,4].astype(bool),fixed_y):
        raise ValueError('boundary flags differ')
    if f['loads'] or f['stages'][:2]!=[[1,'GRAVITY',None,None,1,'自重'],[2,'SRM',None,None,1,'強度低減']]:
        raise ValueError('loading differs')
    return dict(verified=True,sha256=c['source']['sha256'],nodes=661,elements=200,exterior_edges=len(boundary),
                fixed_ux=int(fixed_x.sum()),fixed_uy=int(fixed_y.sum()),accepted=None)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixtures',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.write_text(json.dumps(audit(a.fixtures),indent=2)+'\n','utf-8')
