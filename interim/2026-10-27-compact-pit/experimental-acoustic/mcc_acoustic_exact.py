"""All-direction rank-one criterion for a frozen continuum MCC plastic branch.

This is neither the discrete FEM tangent nor a determination of actual loading.
See exact-criterion.md for the simplex reduction and numerical limitations.
"""
import numpy as np
from mcc_acoustic_diagnostic import acoustic_tensor


def criterion(data, *, relative_decision_margin=1e-10):
    k,g,den=[float(data[n]) for n in ['bulk_kpa','shear_kpa','denominator_kpa3']]
    cn=np.asarray(data['cn'],dtype=float)
    if cn.shape!=(6,) or not np.isfinite(cn).all() or not np.isfinite([k,g,den]).all() or min(k,g,den)<=0:
        raise ValueError('finite positive moduli/consistency denominator required')
    if not np.isfinite(relative_decision_margin) or not 0<relative_decision_margin<=1e-6:
        raise ValueError('bounded positive dimensionless decision margin required')
    a=np.array([[cn[0],cn[3],cn[5]],[cn[3],cn[1],cn[4]],[cn[5],cn[4],cn[2]]])
    principal,axes=np.linalg.eigh(a)
    low,high=float(principal[0]),float(principal[-1]);span=high-low
    eta=g/(k+4*g/3);beta=(k+g/3)/(k+4*g/3)
    mean_opt=float(np.clip((low+high)/(2*beta),low,high))
    fraction=0. if span==0 else float(np.clip((mean_opt-low)/span,0.,1.))
    mean=(1-fraction)*low+fraction*high
    # Variance + eta*mean^2 avoids subtracting nearly equal second moments.
    maximum=fraction*(1-fraction)*span*span+eta*mean*mean
    ratio=float(maximum/(g*den))
    normal=np.sqrt(1-fraction)*axes[:,0]+np.sqrt(fraction)*axes[:,-1]
    normal/=np.linalg.norm(normal)
    elastic=acoustic_tensor(data,normal,plastic=False)
    plastic=acoustic_tensor(data,normal,plastic=True)
    v=a@normal;z=np.linalg.solve(elastic,v)
    independently_contracted_ratio=float(v@z/den)
    if not np.isfinite(ratio) or abs(ratio-independently_contracted_ratio)>2e-10*max(1,ratio):
        raise ValueError('analytic maximum and independent acoustic contraction differ')
    if np.linalg.norm(z)>0:z/=np.linalg.norm(z)
    else:z=axes[:,0]
    rayleigh=float(z@plastic@z)
    eigenvalue=float(np.linalg.eigvalsh(plastic)[0])
    if ratio>1+relative_decision_margin:
        if not rayleigh<0 or not eigenvalue<0:raise ValueError('negative criterion has no verified negative witness')
        decision='negative_continuum_candidate'
    elif ratio<1-relative_decision_margin:decision='positive_continuum_candidate_all_directions'
    else:decision='near_boundary_indeterminate'
    return dict(maximum_rank_one_ratio=ratio,decision=decision,relative_decision_margin=relative_decision_margin,
                extremal_CN_principal_values_kpa2=[low,high],extreme_axis_mixing_fraction=fraction,
                normal=normal.tolist(),polarization=z.tolist(),rayleigh_kpa=rayleigh,
                witness_minimum_eigenvalue_kpa=eigenvalue,
                actual_loading_branch_not_determined=True,is_not_discrete_or_global_tangent=True)
