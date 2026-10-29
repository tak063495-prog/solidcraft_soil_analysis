"""Independent reductions and Jacobian checks before any native comparison."""
from pathlib import Path
import copy
import json
import sys
import unittest
import numpy as np
from reference import Grid, solve, ROOT
from richards_column_reference import solve as column


def case():
    value = json.loads((ROOT/'examples/richards-rain-pause.json').read_text('utf-8'))
    value['analysis']['periods'] = value['analysis']['periods'][:2]
    for i, period in enumerate(value['analysis']['periods']):
        period['end_time_s'] = (i+1)*300.
    return value


class Reference(unittest.TestCase):
    def test_coloured_jacobian_matches_full_dense_finite_differences(self):
        c = case()
        g = Grid(c, 4)
        z, y, x = np.indices((4,4,4))
        h = .12 + .03*np.sin(x+2*y+3*z)
        old = g.field(np.zeros_like(h), 0)[1]
        coeff = g.jacobian(h, old, 50., 3e-6)
        matrix = np.zeros((64,64))
        np.add.at(matrix, (np.arange(64)[:,None],g.columns), coeff)
        dense = np.empty_like(matrix)
        for i in range(64):
            a, b = h.copy(), h.copy()
            delta = 1e-6*(1+abs(h.flat[i]))
            a.flat[i] += delta
            b.flat[i] -= delta
            dense[:,i] = ((g.residual(a,old,50.,3e-6)-g.residual(b,old,50.,3e-6))/(2*delta)).ravel()
        np.testing.assert_allclose(matrix,dense,rtol=1e-8,atol=1e-16)
        self.assertGreater(float(np.max(np.abs(matrix-matrix.T))),1e-8)

    def test_full_width_matches_independently_assembled_1d_and_mass(self):
        c = case()
        actual = solve(c, 4, 50., local=False)
        expected = column(c, 4, 50.)
        for a,b in zip(actual['frames'],expected['frames'],strict=True):
            h = np.asarray(a['head_m'])
            np.testing.assert_allclose(h,np.broadcast_to(np.asarray(b['head_m'])[:,None,None],h.shape),rtol=0,atol=2e-9)
            self.assertLess(abs(a['storage_offset_m3']-b['storage_offset_m']),1e-10)
            self.assertLess(a['transverse_flux_rms_m_s'],1e-14)
        self.assertLess(actual['maximum_integrated_mass_error_m3'],1e-10)

    def test_local_patch_has_transverse_flow_and_closed_global_balance(self):
        result = solve(case(), 4, 50.)
        self.assertEqual(result['rain_patch_area_m2'],.25)
        self.assertGreater(result['frames'][0]['transverse_flux_rms_m_s'],1e-12)
        self.assertGreater(np.ptp(np.asarray(result['frames'][0]['head_m'])[-1]),1e-5)
        self.assertLess(result['maximum_integrated_mass_error_m3'],1e-10)
        for frame in result['frames']:
            self.assertGreaterEqual(frame['cumulative_surface_routed_water_m3'],-1e-14)
        self.assertEqual(result['frames'][0]['cumulative_precipitation_m3'],.25*3e-6*300)


if __name__ == '__main__':
    unittest.main()
