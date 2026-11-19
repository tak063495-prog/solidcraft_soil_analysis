"""One original 12-Tet4 Newton step, with binary64 NumPy and portable state."""
from decimal import localcontext
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

HERE = Path(__file__).resolve().parents[1]/'interim/2027-01-12-a04-fixed-order-portable-reference'


class PortableHydraulics(unittest.TestCase):
    def test_first_vg_wetting_step(self):
        path = sys.path[:]
        try:
            sys.path.insert(0,str(HERE))
            from a04ref_matrix import driver
            from a04ref_decimal import PRECISION, model_type
            bound = driver.bind()
        finally:
            sys.path[:] = path
        values = bound['study'].fixture(4,False)
        log = []
        with patch.object(np,'longdouble',np.float64), localcontext() as ctx:
            ctx.prec = PRECISION
            self.assertEqual(np.finfo(np.longdouble).nmant,52)
            m = model_type(bound['soil_reference'].Model,log)(*values)
            v = m.initial_hydraulic()
            p = m.flow['analysis']['periods'][1]
            result = m.full.step(v['h'],v['trace'],v['storage'],p,.25)
        self.assertGreater(result['iterations'],0)
        self.assertLess(result['maximum_scaled_residual'],.01)
        self.assertLessEqual(result['maximum_true_linear_residual'],1e-10)
        self.assertEqual(len(log),result['iterations']+1)
        for row in log[:-1]:
            self.assertTrue(row['trials'])
            self.assertLess(float(row['trials'][-1]['norm']),float(row['norm']))
        self.assertEqual(result['outflow_m3_s'].shape,(12,4))
