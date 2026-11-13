"""Trace/load algebra only, not a JOINT constitutive or BIRTH implementation."""
import numpy as np
from reference import Darcy


def trace(kind, z, normal, head_a=None, head_b=None, level=None, soil_side=None, gamma_w=9.81):
    n = np.asarray(normal, dtype=float)
    if n.shape != (2,) or not np.isclose(np.linalg.norm(n), 1., rtol=0., atol=1e-14):
        raise ValueError('explicit unit a-to-b normal required')
    if kind == 'continuous':
        if head_a is None or head_b is None or head_a != head_b:
            raise ValueError('CONTINUOUS requires explicit identical hydraulic trace; no averaging')
        p = gamma_w*(head_a-z)
        pore_a = pore_b = p
    elif kind == 'soil':
        if soil_side not in ('a', 'b'):
            raise ValueError('SOIL side required')
        h = head_a if soil_side == 'a' else head_b
        if h is None:
            raise ValueError('missing SOIL head')
        p = gamma_w*(h-z)
        pore_a, pore_b = (p, 0.) if soil_side == 'a' else (0., p)
    elif kind == 'prescribed':
        if level is None or any(h is not None and h != level for h in (head_a, head_b)):
            raise ValueError('explicit prescribed head required; conflicting heads rejected')
        p = gamma_w*(level-z)
        pore_a = pore_b = p
    else:
        raise ValueError('unknown policy')
    if not np.isfinite(p) or p < 0:
        raise ValueError('undefined/negative saturated pressure; no clipping')
    return dict(pressure_kpa=p, pore_a_kpa=pore_a, pore_b_kpa=pore_b,
                traction_a=(-p*n).tolist(), traction_b=(p*n).tolist())


def reservoir_case():
    """PRESCRIBED internal reservoir: its water exchange must not disappear.

    Different K on both sides and distinct end heads; each subdomain is source
    free. The reservoir is external to both subdomains, unlike CONTINUOUS.
    """
    low = Darcy(edges=(0., 1.2), conductivity=(1e-5,), bottom_head=4.4, top_head=4.)
    high = Darcy(edges=(1.2, 3.), conductivity=(4e-5,), bottom_head=4., top_head=3.6)
    bottom = -low.flux
    top = high.flux
    reservoir_out = low.flux-high.flux
    return dict(low=low, high=high, bottom_out=bottom, top_out=top,
                reservoir_out=reservoir_out, water_defect=bottom+top+reservoir_out)


def mcc_invariants(state):
    """Only documented material-point definitions, not an initialization rule."""
    stress = np.asarray(state['effective_stress_kpa'], dtype=float)
    if stress.shape != (6,) or not np.all(np.isfinite(stress)):
        raise ValueError('explicit six-component MCC effective stress required')
    pc, e, m = (float(state[key]) for key in ('pc_kpa', 'void_ratio', 'M'))
    p = float(sum(stress[:3])/3)
    if not np.all(np.isfinite([pc, e, m])) or min(p, pc, e, m) <= 0:
        raise ValueError('MCC requires positive p/pc/e/M; no stress-free fallback')
    deviator = stress[:3]-p
    q2 = float(1.5*(np.dot(deviator, deviator)+2*np.dot(stress[3:], stress[3:])))
    return dict(p_kpa=p, q_kpa=float(np.sqrt(q2)), yield_function_kpa2=q2+m*m*p*(p-pc),
                initialization='仕様未定義', accepted=None)
