"""Fixed-order reductions for the C-02 observation, never acceptance corrections."""
import math
import numpy as np
from c02diag_transfer import norm,dot


def difference(a,b,weights=None):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    if a.shape!=b.shape:
        raise ValueError('history shape mismatch')
    absolute=norm(a-b,weights);reference=norm(b,weights)
    return dict(absolute_l2=absolute,reference_l2=reference,relative_l2=absolute/max(reference,1e-12),maximum_absolute=float(np.max(abs(a-b))))


def triple(values,weights=None):
    a,b,c=values;d0=a-b;d1=b-c;n0=norm(d0,weights);n1=norm(d1,weights)
    return dict(coarse_middle=difference(a,b,weights),middle_fine=difference(b,c,weights),order=math.log2(n0/n1) if min(n0,n1)>0 else None,difference_cosine=dot((d0 if weights is None else d0*np.broadcast_to(weights,d0.shape)).ravel(order='C'),d1.ravel(order='C'))/(n0*n1) if min(n0,n1)>0 else None)


def time_identity(a0,b0,a1,b1,labels):
    """S(t1)-S(t0)=T(h0)-T(h1), aligned physical points and accepted times."""
    s0=a0-b0;s1=a1-b1;change=s1-s0;rhs=(a1-a0)-(b1-b0)
    label=np.asarray(labels,dtype=bool)
    masks={'all_closed':(~label).all(axis=0),'all_open':label.all(axis=0),'label_disagreement':label.any(axis=0)&~label.all(axis=0)}
    total=norm(change)**2
    split={}
    for name,mask in masks.items():
        squared=norm(change[mask])**2
        split[name]=dict(point_times=int(mask.sum()),change_absolute_l2=math.sqrt(squared),fraction_change_squared=squared/total if total else None,fine_spatial_absolute_l2=norm(s1[mask]))
    return dict(coarse_spatial_absolute_l2=norm(s0),fine_spatial_absolute_l2=norm(s1),spatial_change_absolute_l2=norm(change),change_over_fine_spatial=norm(change)/max(norm(s1),1e-12),identity_roundoff_absolute_l2=norm(change-rhs),label_partition=split)
