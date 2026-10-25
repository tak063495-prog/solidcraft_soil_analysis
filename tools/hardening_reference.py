"""Independent IDH checks, analytic shear elimination and scalar cap brackets.

No Rust tangents, trial stresses or hardening histories generate the reference.
Eliminate p,q using eta=lambda_cap/radius; shear consistency is a quadratic.
Only the cap needs a scalar bracket. Rust enumerates 4-variable Newton sets.
"""
import math
import numpy as np


def invariants(stress):
    stress=np.asarray(stress);p=stress[..., :3].mean(axis=-1)
    dev=stress.copy();dev[..., :3]-=p[...,None]
    q=np.sqrt(1.5*(dev[...,:3]**2).sum(axis=-1)+3*(dev[...,3:]**2).sum(axis=-1))
    return p,q


def shear_ratio(par,kappa):
    return par['shear_ratio_initial']+(par['shear_ratio_ultimate']-par['shear_ratio_initial'])*kappa/(par['shear_reference_plastic_strain']+kappa)


def cap_pressure(par,initial,volume):
    power=1-par['stiffness_exponent'];b=par['cap_hardening_modulus_ref_kpa']/par['reference_pressure_kpa']
    if power==0:return initial*np.exp(b*volume)
    return initial*np.exp(np.log1p(power*b*volume/(initial/par['reference_pressure_kpa'])**power)/power)


def audit(par,states):
    stress=np.asarray([s['effective_stress_compression_positive_kpa'] for s in states]);p,q=invariants(stress)
    shear=np.asarray([s['shear_plastic_strain'] for s in states]);volume=np.asarray([s['cap_plastic_volumetric_strain'] for s in states])
    pc=np.asarray([s['cap_pressure_kpa'] for s in states]);initial=np.asarray([s['initial_cap_pressure_kpa'] for s in states])
    expected=cap_pressure(par,initial,volume)
    assert np.all(p>par['minimum_mean_effective_pressure_kpa']) and np.all(shear>=0) and np.all(volume>=0)
    assert np.allclose(pc,expected,rtol=5e-9,atol=1e-9)
    trace=np.asarray([sum(s['plastic_strain_compression_positive'][:3]) for s in states])
    assert np.allclose(trace,volume-par['dilation_ratio']*shear,rtol=1e-9,atol=1e-9)
    assert all(s['plastic_work_kpa']>=0 and s['void_ratio']>0 for s in states)
    fs=(q-shear_ratio(par,shear)*(p+par['shear_pressure_shift_kpa']))/np.maximum(par['reference_pressure_kpa'],par['shear_ratio_ultimate']*(p+par['shear_pressure_shift_kpa']))
    fc=(np.hypot(p,q/par['cap_axis_ratio'])-pc)/np.maximum(pc,par['reference_pressure_kpa'])
    return np.maximum(fs,fc)


def substep(par,old,de):
    stress=np.asarray(old['effective_stress_compression_positive_kpa']);de=np.asarray(de)
    po,qo=invariants(stress);e=par['unloading_young_ref_kpa']*(old['cap_pressure_kpa']/par['reference_pressure_kpa'])**par['stiffness_exponent']
    bulk=e/(3*(1-2*par['poisson']));g=e/(2*(1+par['poisson']))
    trial=stress+np.r_[2*g*de[:3]+(bulk-2*g/3)*de[:3].sum(),g*de[3:]]
    pt,qt=invariants(trial)
    axis=par['cap_axis_ratio'];dilation=par['dilation_ratio'];oldk=old['shear_plastic_strain']
    # For fixed eta, p=A+B*ls, q=C-D*ls and
    # M=Mmax-(Mmax-M0)*gamma_ref/(gamma_ref+k_old+ls).
    # Multiplication by the positive hardening denominator gives a quadratic
    # with exactly one positive root when the elastic trial violates shear.
    def state_at(eta):
        def values(ls):
            p=(pt+bulk*dilation*ls)/(1+bulk*eta)
            q=(qt-3*g*ls)/(1+3*g*eta/axis**2)
            return p,q,q-shear_ratio(par,oldk+ls)*(p+par['shear_pressure_shift_kpa'])
        p,q,fs=values(0.);ls=0.
        if fs>1e-12*max(1.,pt,qt):
            aa=pt/(1+bulk*eta);bb=bulk*dilation/(1+bulk*eta)
            cc=qt/(1+3*g*eta/axis**2);dd=3*g/(1+3*g*eta/axis**2)
            h=par['shear_reference_plastic_strain']+oldk
            z=(par['shear_ratio_ultimate']-par['shear_ratio_initial'])*par['shear_reference_plastic_strain']
            c0=cc-par['shear_ratio_ultimate']*(aa+par['shear_pressure_shift_kpa'])
            c1=-dd-par['shear_ratio_ultimate']*bb
            a=c1;b=c0+c1*h+z*bb;c=c0*h+z*(aa+par['shear_pressure_shift_kpa'])
            discriminant=math.sqrt(b*b-4*a*c)
            ls=(-b-discriminant)/(2*a) if b>=0 else 2*c/(-b+discriminant)
            assert 0<=ls<=qt/(3*g)*(1+1e-12)
            p,q,fs=values(ls)
        radius=math.hypot(p,q/axis);volume=eta*p
        pc=cap_pressure(par,old['cap_pressure_kpa'],volume)
        return p,q,ls,radius,volume,float(pc),radius-pc
    eta=0.;p,q,ls,radius,volume,pc,fc=state_at(0.)
    if fc>1e-12*max(pt,qt,old['cap_pressure_kpa'],1.):
        lo=0.;hi=1/max(bulk,g)
        for _ in range(80):
            if state_at(hi)[-1]<0:break
            hi*=2
        else:raise AssertionError('no independent positive cap bracket')
        for _ in range(58):
            mid=(lo+hi)/2
            if state_at(mid)[-1]>0:lo=mid
            else:hi=mid
        eta=(lo+hi)/2;p,q,ls,radius,volume,pc,fc=state_at(eta)
    assert p>par['minimum_mean_effective_pressure_kpa']
    dev=trial-np.r_[np.repeat(pt,3),np.zeros(3)]
    nq=dev*np.array([1.5,1.5,1.5,3,3,3])/qt if qt>0 else np.zeros(6)
    flow=ls*(nq-np.r_[np.repeat(dilation/3,3),np.zeros(3)])+eta*(np.r_[np.repeat(p/3,3),np.zeros(3)]+q*nq/axis**2)
    result=np.r_[np.repeat(p,3),np.zeros(3)]+(dev*q/qt if qt>0 else np.zeros(6))
    return {'effective_stress_compression_positive_kpa':result.tolist(),
        'plastic_strain_compression_positive':(np.asarray(old['plastic_strain_compression_positive'])+flow).tolist(),
        'void_ratio':(1+old['void_ratio'])*math.exp(-de[:3].sum())-1,
        'shear_plastic_strain':oldk+ls,'cap_plastic_volumetric_strain':old['cap_plastic_volumetric_strain']+volume,
        'initial_cap_pressure_kpa':old['initial_cap_pressure_kpa'],'cap_pressure_kpa':pc,
        'plastic_work_kpa':old['plastic_work_kpa']+max(0.,float(result@flow)),
        'last_return':'shear_cap' if ls>0 and eta>0 else ('shear' if ls>0 else ('cap' if eta>0 else 'elastic'))}


def integrate(par,old,strain,count):
    state=old
    for _ in range(count):state=substep(par,state,np.asarray(strain)/count)
    assert max(audit(par,[state]))<2e-8
    return state
