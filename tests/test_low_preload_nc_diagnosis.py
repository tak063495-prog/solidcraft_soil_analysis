"""Frozen-checkpoint diagnostics and refusal of changed release conditions."""
import copy
import json
import os
from pathlib import Path
import sys
import unittest
import zipfile

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
import probe_low_preload_nc as low
from diagnose_low_preload_nc import affine_load_check, state_diagnostics
from low_preload_nc_reference import ClosedStickSystem
from verify_high_order import physical


class LowNCDiagnosis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh, cls.ground, cls.case = low.expected_inputs(1/128, minimum_increment=1/16384, dense_lu=True)
        archive = ROOT/'interim/2026-11-08-low-preload-nc-followup/native-evidence.zip'
        with zipfile.ZipFile(archive) as z:
            cls.saved_run = json.loads(z.read('dense-1-128-min16384-v2/job/recovery-checkpoint-1.json'))['run']

    def test_saved_fields_exclude_live_contact_switch_and_retain_history(self):
        before = copy.deepcopy(self.saved_run)
        report = state_diagnostics(self.mesh, self.ground, self.case, self.saved_run)
        self.assertEqual(report['active_contact_points'], 0)
        self.assertEqual(report['active_material_points'], 162)
        self.assertEqual(len(report['continuum_singular_values_kn_m']), 54)
        self.assertTrue(report['history_and_force']['removed_joint_and_material_history_frozen'])
        self.assertLess(report['history_and_force']['last_committed_force_check']['relative_force_residual'], 1e-8)
        self.assertGreater(report['affine_release']['remaining_holding_norm_kn'], 0)
        self.assertEqual(self.saved_run, before)

    def test_changed_free_load_is_refused_instead_of_being_named_equilibrium(self):
        state = self.saved_run['last_committed']
        model = ClosedStickSystem(self.mesh, self.ground, state['active_elements'], self.case)
        holding = np.asarray(self.saved_run['excavation_audits'][0]['initial_holding_load_kn'])
        endload = physical(self.mesh, self.ground, np.asarray(state['active_elements']), 0., True,
                           joint_system=model.joint)[0].ravel()
        changed = copy.deepcopy(state)
        # A free nodal load, not a constrained reaction or inactive node.
        changed['applied_nodal_load_kn'][model.free[0]] += 1.
        with self.assertRaisesRegex(ValueError, 'affine release'):
            affine_load_check(model, changed, endload, holding)

    def test_difference_width_refines_only_existing_numerical_control(self):
        original = low.expected_inputs()
        changed = low.expected_inputs(tangent_perturbation=1e-9)
        changed[2]['settings']['tangent_perturbation'] = 1e-7
        self.assertEqual(changed, original)
        for width in (0., 1e-10, 1e-6, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                low.expected_inputs(tangent_perturbation=width)


if __name__ == '__main__':
    unittest.main()
