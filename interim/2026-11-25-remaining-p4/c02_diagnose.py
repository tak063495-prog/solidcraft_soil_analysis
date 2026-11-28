"""Matched-degree and matched-grid comparisons on the 448 P3 anchors."""
from pathlib import Path
import sys
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
out=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'c02'
P2=HERE.parent/'2026-11-25-remaining5-p2/c02'
P3=HERE.parent/'2026-11-25-remaining-p3/c02'
paths={'ref':P2/'p3-d8-n256','ref4':P2/'p4-d8-n256','ref512':P2/'p3-d8-n512',
       '3a':P3/'p3','4b':P3/'p4','5c':P3/'p5'}
for label,tag in [('4a','p4n256'),('5a','p5n256'),('5b','p5n512'),('3b','p3n512')]:
    paths[label]=out/tag
data={k:np.load(v/'principal.npz') for k,v in paths.items()}
ref=data['ref']
w=ref['weight']
xy=ref['position'][:,:2]
distance=np.minimum(xy,1-xy).min(axis=1)
for d in data.values():
    for key in ('parent','weight','position'):
        np.testing.assert_allclose(ref[key],d[key],atol=2e-15,rtol=0)


def compare(a,b):
    x,y=data[a],data[b]
    row={}
    for field in ('slip','signed_area_integral'):
        scale=np.sqrt(w)[None,:,None] if field=='slip' else 1.
        diff=(x[field]-y[field])*scale
        base=ref[field]*scale
        sq=diff**2
        row[field]=dict(percent_ref_norm=100*np.linalg.norm(diff)/np.linalg.norm(base),
                        stage_squared_shares=[np.sum(v)/np.sum(sq) for v in np.split(sq,7)])
        if field=='slip':
            row[field]['edge_squared_shares']={str(t):sq[:,distance<t].sum()/sq.sum() for t in (.05,.1,.2)}
    return row


pairs=[('3a','ref'),('4b','ref'),('5c','ref'),('4a','3a'),('5a','4a'),('4b','3b'),('5b','4b'),
       ('3b','3a'),('4b','4a'),('5b','5a'),('5c','5b'),('4a','ref4'),('ref4','ref'),('ref512','ref')]
report=dict(denominator='fixed P3 d8 N256 reference for each field',edge_area={str(t):w[distance<t].sum()/w.sum() for t in (.05,.1,.2)},
            inputs={k:dict(path=str(v.relative_to(HERE.parents[1]) if v.is_relative_to(HERE.parents[1]) else v),sha256=hashlib.sha256((v/'principal.npz').read_bytes()).hexdigest()) for k,v in paths.items()},
            pairs={f'{a}-{b}':compare(a,b) for a,b in pairs})
(out/'comparisons.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
