"""Continuum PDE control for C-04, independent of RT0 algebra agreement."""
import copy
import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('c04_series_tests', ROOT/'interim/2026-11-22-wet-joint-transient-water/saturated_control.py')
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)


def fixture():
    case = s.s.hydraulic_case(4); case['initial_head'] = dict(head_m=2.)
    for p in case['analysis']['periods']:
        value = 2. if p['end_time_s'] <= 1. else 2.2
        for b in p['boundaries']: b['condition']['value']['head_m'] = value
    return case


class TransientSeries(unittest.TestCase):
    def test_series_satisfies_original_diffusion_equation_and_instantaneous_flux(self):
        ref = s.Series(fixture()); z = np.linspace(.1, .9, 7); t = 1.2; h = 1e-4
        H, q, _ = ref.fields(z, t)
        dt = (ref.fields(z, t+h)[0]-ref.fields(z, t-h)[0])/(2*h)
        dz2 = (ref.fields(z+h, t)[0]-2*H+ref.fields(z-h, t)[0])/h**2
        dq = -ref.k*(ref.fields(z+h, t)[0]-ref.fields(z-h, t)[0])/(2*h)
        np.testing.assert_allclose(dt, ref.diffusivity*dz2, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(q, dq, rtol=1e-6, atol=1e-9)
        self.assertAlmostEqual(ref.fields(np.array([0., 1.]), 1.)[0][0], 2.)
        np.testing.assert_array_equal(ref.fields(np.array([0., 1.]), 1.25)[0], [2.2, 2.2])

    def test_series_storage_rate_equals_boundary_inflow_and_refuses_unsaturated_substitution(self):
        case = fixture(); ref = s.Series(case); t = 1.2; h = 1e-5
        rate = (ref.fields(np.array([.5]), t+h)[2]-ref.fields(np.array([.5]), t-h)[2])/(2*h)
        flux = ref.fields(np.array([0., 1.]), t)[1]
        self.assertAlmostEqual(rate, flux[0]-flux[1], places=10)
        bad = copy.deepcopy(case); bad['initial_head']['head_m'] = .3
        with self.assertRaisesRegex(ValueError, 'strict saturation'): s.Series(bad)


if __name__ == '__main__': unittest.main()
