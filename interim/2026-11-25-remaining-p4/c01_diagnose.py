"""Factorial route/budget comparison and contact-transition counts at common anchors."""
from pathlib import Path
import json
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
P2=HERE.parent/'2026-11-25-remaining5-p2'
P3=HERE.parent/'2026-11-25-remaining-p3'
sys.path.insert(0,str(P2))
from c01_score import observed,STAGES,CATEGORIES
work=Path(sys.argv[1])
out=Path(sys.argv[2]) if len(sys.argv)>2 else HERE/'c01'
out.mkdir()
sizes=(256,512,1024)


def field(folder,name):
    meta=json.loads((folder/'fields/manifest.json').read_bytes())['fields'][name]
    return np.fromfile(folder/'fields'/(name+'.bin'),dtype=meta['dtype']).reshape(1792,*meta['shape'])


def events(folder):
    s=field(folder,'contact_status').reshape(7,256,-1)
    rows={}
    for i,(stage,values) in enumerate(zip(STAGES,s)):
        before=np.concatenate(([s[i-1,-1] if i else values[0]],values[:-1]))
        pairs={}
        for a in range(3):
            for b in range(3):
                if a==b:
                    continue
                t,point=np.nonzero((before==a)&(values==b))
                if len(t):
                    pairs[f'{a}->{b}']=dict(count=len(t),first_anchor=int(t.min()+1),last_anchor=int(t.max()+1))
        rows[stage]=pairs
    return rows


original=[work/f'original{n}' for n in sizes]
alt=work/'alternate256'
alternate_manifest=json.loads((alt/'fields/manifest.json').read_bytes())
assert alternate_manifest==json.loads((P3/'c01/fixed/manifest-256.json').read_bytes())
names=json.loads((original[0]/'fields/manifest.json').read_bytes())['fields']
fields={}
for name in names:
    a,b,c=[field(path,name) for path in original]
    fields[name]=observed(a,b,c) if name not in CATEGORIES else dict(mismatch=[int((a!=b).sum()),int((b!=c).sum())])
slip=[field(path,'plastic_slip') for path in original]
archives={k:np.load(v/'slip-common.npz') for k,v in [('orig_scaled',P2/'c01'),('alt_fixed',P3/'c01/fixed'),('alt_scaled',P3/'c01/scaled')]}
np.testing.assert_array_equal(slip[-1],archives['orig_scaled']['n1024'])
np.testing.assert_array_equal(field(alt,'plastic_slip'),archives['alt_fixed']['n256'])
report=dict(original_fixed_fields=fields,events_256={'original':events(original[0]),'alternate':events(alt)},
            stages=list(STAGES),routes={})
for label,series in [('orig_fixed',slip),*((key,[f[f'n{n}'] for n in sizes]) for key,f in archives.items())]:
    summary={}
    if label=='orig_fixed':
        routes={str(n):json.loads((d/'result.json').read_bytes()) for n,d in zip(sizes,original)}
    else:
        folder={'orig_scaled':P2/'c01','alt_fixed':P3/'c01/fixed','alt_scaled':P3/'c01/scaled'}[label]
        routes=json.loads((folder/'summary.json').read_bytes())['routes']
    for stage_id,stage in enumerate(STAGES):
        idx=slice(stage_id*256,(stage_id+1)*256)
        a,b,c=[x[idx] for x in series]
        errors=[a-b,b-c]
        summary[stage]=dict(**observed(a,b,c),new_error_norm=[float(np.linalg.norm(x-x[0])) for x in errors],
                            relative_to_recontact_endpoint=[float(np.linalg.norm(x-(series[j][1023]-series[j+1][1023]))/np.linalg.norm(x)) if np.linalg.norm(x)>0 else None for j,x in enumerate(errors)],
                            attempts=[routes[str(n)]['stages'][stage_id]['attempts'] for n in sizes],
                            rebalances=[routes[str(n)]['stages'][stage_id]['reequilibrated_coarse_mid_fine'] for n in sizes])
    report['routes'][label]=dict(overall=observed(*series),by_stage=summary)
# Resolve signed reversals against the shared shear endpoint, at the existing slip budget.
axis=slip[0][511].sum(axis=0)
axis/=np.linalg.norm(axis)
report['shear_axis']=axis.tolist()
report['signed_activity_256']={}
for name,history in [('original',slip[0]),('alternate',archives['alt_fixed']['n256'])]:
    change=np.diff(history,axis=0,prepend=np.zeros_like(history[:1]))@axis
    per_stage={}
    for stage,chunk in zip(STAGES,np.split(change,7)):
        entry={}
        for label,sign in [('positive',1),('negative',-1)]:
            anchors,points=np.where(sign*chunk>(1e-4/1024**2))
            if len(anchors):
                entry[label]=dict(point_anchors=len(anchors),first=int(anchors.min()+1),last=int(anchors.max()+1))
        per_stage[stage]=entry
    report['signed_activity_256'][name]=per_stage
common_prefix={str(n):np.max(abs(x[:1024]-archives['alt_fixed'][f'n{n}'][:1024])) for n,x in zip(sizes,slip)}
report['original_alternate_fixed_first_four_stage_max_difference_m']=common_prefix
for path in original+[alt]:
    (out/f'{path.name}-manifest.json').write_bytes((path/'fields/manifest.json').read_bytes())
    for filename in ('input.json','result.json'):
        (out/f'{path.name}-{filename}').write_bytes((path/filename).read_bytes())
np.savez_compressed(out/'original-fixed-slip.npz',**{f'n{n}':x for n,x in zip(sizes,slip)})
(out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
