"""Deterministic observation metrics. Original L2 definition, no new pass gate."""
import math
from pathlib import Path
import json
import numpy as np

CATEGORIES=('material_return','contact_status','open')


def norm(values):
    return math.sqrt(math.fsum(float(x)*float(x) for x in np.asarray(values).ravel(order='C')))


def dot(a,b):
    a=np.asarray(a).ravel(order='C');b=np.asarray(b).ravel(order='C')
    if a.shape!=b.shape:raise ValueError('vector correspondence differs')
    return math.fsum(float(x)*float(y) for x,y in zip(a,b))


def cosine(a,b):
    den=norm(a)*norm(b)
    return dot(a,b)/den if den else None


def load(path):
    rows=[json.loads(s) for s in Path(path).read_text(encoding='utf-8').splitlines()]
    keys=sorted(rows[0]['fields'])
    if any(sorted(row['fields'])!=keys for row in rows):raise ValueError('field inventory differs')
    result={k:np.asarray([row['fields'][k] for row in rows]) for k in keys}
    if not all(np.isfinite(v).all() for v in result.values()):raise ValueError('nonfinite field')
    return [row['row'] for row in rows],result


def metric(a,b):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('invalid comparison')
    d=a-b;absolute=norm(d);reference=norm(b);ix=np.unravel_index(int(np.argmax(np.abs(d))),d.shape)
    return dict(absolute_l2=absolute,reference_l2=reference,relative_l2=absolute/max(reference,1e-12),
                maximum_absolute=float(np.max(np.abs(d))),maximum_index=list(map(int,ix)),
                signed_sum=math.fsum(map(float,d.ravel(order='C'))),components=d.size)


def compare(a,b):
    ra,fa=a;rb,fb=b
    if [(r['stage'],r['progress']) for r in ra]!=[(r['stage'],r['progress']) for r in rb]:raise ValueError('history alignment differs')
    if set(fa)!=set(fb):raise ValueError('field inventory differs')
    fields={};phases={};endpoints={};indices=[i for i,r in enumerate(ra) if r['progress']==1.]
    for k in sorted(fa):
        if k in CATEGORIES:
            fields[k]=dict(mismatches=int((fa[k]!=fb[k]).sum()),components=fa[k].size)
        else:
            fields[k]=metric(fa[k],fb[k]);endpoints[k]=metric(fa[k][indices],fb[k][indices])
            for stage in dict.fromkeys(r['stage'] for r in ra):
                ix=[i for i,r in enumerate(ra) if r['stage']==stage]
                phases.setdefault(stage,{})[k]=metric(fa[k][ix],fb[k][ix])
    return dict(accepted_pairs=len(ra),all_fields=fields,endpoints=endpoints,stages=phases,
                numerical_acceptance=False)


def budget(original,contrast):
    """Vector change is descriptive, never an additive causal percentage."""
    return dict(original_absolute_l2=norm(original),contrast_absolute_l2=norm(contrast),
                contrast_over_original=norm(contrast)/max(norm(original),1e-30),
                difference_over_original=norm(contrast-original)/max(norm(original),1e-30),
                cosine=cosine(original,contrast),
                projection_on_original=dot(contrast,original)/max(dot(original,original),1e-60))
