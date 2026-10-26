"""Reconstruct compact MCC seed histories from original coefficients/geometry."""
import numpy as np
from mcc_rate_fem import seed as independent_seed


def verify_affine_seed(mesh,ground,states):
    if ground['initialization']['geostatic']['seed']['kind']!='mcc_affine_fields':
        raise ValueError('affine MCC definition required')
    expected=independent_seed(mesh,ground).reshape(-1,14)
    keys={'effective_stress_compression_positive_kpa','preconsolidation_pressure_kpa',
          'void_ratio','plastic_strain_compression_positive'}
    if len(states)!=len(expected) or any(set(s)!=keys for s in states):
        raise ValueError('affine MCC seed history inventory/fields differ')
    for name,values,absolute in (
        ('effective_stress_compression_positive_kpa',expected[:,:6],1e-10),
        ('preconsolidation_pressure_kpa',np.exp(expected[:,6]),1e-10),
        ('void_ratio',np.expm1(expected[:,7]),1e-12)):
        actual=np.asarray([s[name] for s in states],dtype=float)
        if (actual.shape!=values.shape or not np.isfinite(actual).all()
                or not np.allclose(actual,values,rtol=1e-12,atol=absolute)):
            raise ValueError('affine MCC initial history differs from input coefficients/NCL: '+name)
    plastic=np.asarray([s['plastic_strain_compression_positive'] for s in states],dtype=float)
    if plastic.shape!=expected[:,8:].shape or not np.array_equal(plastic,expected[:,8:]):
        raise ValueError('affine MCC seed plastic history is not zero')
