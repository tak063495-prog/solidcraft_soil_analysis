"""Opening-time labels on saved linear gaps; telescoping accounting, not causality."""
import math
import numpy as np
from c02diag_transfer import norm,dot


def catalog(rows,gaps):
    gaps=np.asarray(gaps,dtype=float);result={};previous=np.zeros(gaps.shape[1])
    for stage in dict.fromkeys(r['stage'] for r in rows):
        ids=[i for i,r in enumerate(rows) if r['stage']==stage]
        intervals=[[] for _ in previous];before=previous.copy();t0=0.
        for i in ids:
            t1=float(rows[i]['progress']);after=gaps[i]
            for gp,(g0,g1) in enumerate(zip(before,after,strict=True)):
                opened=None
                if g0>0 and g1>0:opened=(t0,t1)
                elif (g0>0)!=(g1>0):
                    root=t0+(t1-t0)*float(g0)/(float(g0)-float(g1))
                    opened=(t0,root) if g0>0 else (root,t1)
                if opened is not None and opened[0]<opened[1]:
                    old=intervals[gp]
                    if old and old[-1][1]==opened[0]:old[-1]=(old[-1][0],opened[1])
                    else:old.append(opened)
            before=after;t0=t1
        result[stage]=intervals;previous=before
    return result


def clipped(intervals,lo,hi):
    return [(max(a,lo),min(b,hi)) for a,b in intervals if a<hi and b>lo]


def masks(catalogs,stages,steps):
    count=len(catalogs[0][stages[0]]);labels=[];durations=[]
    for stage in stages:
        for j in range(1,steps+1):
            lo=(j-1)/steps;hi=j/steps;row=[];duration=[]
            for gp in range(count):
                pieces=[clipped(c[stage][gp],lo,hi) for c in catalogs]
                row.append(any(x!=pieces[0] for x in pieces[1:]))
                # Pairwise symmetric open-set length, largest pair for 4 paths.
                largest=0.
                for a in range(len(pieces)):
                    for b in range(a+1,len(pieces)):
                        x,y=pieces[a],pieces[b]
                        overlap=math.fsum(max(0.,min(d,f)-max(c,e)) for c,d in x for e,f in y)
                        length=math.fsum(d-c for c,d in x)+math.fsum(d-c for c,d in y)-2*overlap
                        largest=max(largest,length)
                duration.append(largest)
            labels.append(row);durations.append(duration)
    return np.asarray(labels),np.asarray(durations)


def budget(delta,mask,stages,steps):
    delta=np.asarray(delta,dtype=float)
    if delta.shape[:2]!=mask.shape:raise ValueError('opening label shape mismatch')
    previous=np.zeros_like(delta[0]);acc=np.zeros_like(previous);labels=[]
    for value,row in zip(delta,mask,strict=True):
        increment=value-previous;shaped=row.reshape((len(row),)+(1,)*(increment.ndim-1))
        acc=acc+np.where(shaped,increment,0.)
        labels.append(acc.copy());previous=value
    label=np.stack(labels);other=delta-label;total=dot(delta.ravel(),delta.ravel())
    def inner(a,b):return dot(a.ravel(order='C'),b.ravel(order='C'))
    result=dict(total_absolute_l2=math.sqrt(total),opening_label_point_times=int(mask.sum()),opening_label_absolute_l2=norm(label),other_path_absolute_l2=norm(other),opening_label_signed_projection=inner(label,delta)/total if total else None,other_path_signed_projection=inner(other,delta)/total if total else None,cross_inner_product=inner(label,other),identity_residual_l2=norm(delta-(label+other)),stages={})
    for j,s in enumerate(stages):
        d=delta[j*steps:(j+1)*steps];l=label[j*steps:(j+1)*steps];o=other[j*steps:(j+1)*steps]
        ss=inner(d,d)
        result['stages'][s]=dict(difference_squared=ss,opening_label_signed_projection=inner(l,d)/ss if ss else None,other_path_signed_projection=inner(o,d)/ss if ss else None)
    return result,dict(opening_label=label,other_path=other)


def direction_budget(left,right):
    scale=norm(left['opening_label']+left['other_path'])*norm(right['opening_label']+right['other_path'])
    return {a+'/'+b:dot(left[a].ravel(),right[b].ravel())/scale if scale else None for a in left for b in right}
