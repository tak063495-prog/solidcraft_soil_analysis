"""Refuse changed integration, preload and binary before preparing a restart."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import probe_nc_midpoint as nc
from probe_low_preload_nc import sha, write


class MidpointResumeGuards(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name); self.job = self.folder/'old'
        (self.job/'inputs').mkdir(parents=True); (self.job/'solver').mkdir()
        # This fixture tests driver preflight only. The real CLI, not this
        # helper, validates checkpoint/job seals and executable provenance.
        self.binary = self.job/'solver/soil-fem-3d.exe'; self.binary.write_bytes(b'fixture')
        self.digest = sha(self.binary)
        write(self.job/'job.json', {'binding': {'identity': {'solver_binary_sha256': self.digest}}})
        for name, value in zip(('mesh', 'ground', 'case'), nc.inputs('midpoint')):
            write(self.job/'inputs'/(name+'.json'), value)

    def test_changed_scheme_and_100_kpa_preload_are_rejected_before_output(self):
        nc.resume_inputs(self.job, 'midpoint', self.digest)
        with self.assertRaisesRegex(ValueError, 'input mismatch: stages'):
            nc.run(self.binary, self.folder/'wrong-scheme', 'backward_euler', 1., self.job)
        self.assertFalse((self.folder/'wrong-scheme').exists())
        path = self.job/'inputs/case.json'; value = json.loads(path.read_text())
        value['boundary']['surface_loads'][0]['traction_kpa'][2] = -100.
        path.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, 'input mismatch: stages'):
            nc.run(self.binary, self.folder/'wrong-preload', 'midpoint', 1., self.job)
        self.assertFalse((self.folder/'wrong-preload').exists())

    def test_changed_binary_and_binding_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'identical sealed solver'):
            nc.resume_inputs(self.job, 'midpoint', '0'*64)
        self.binary.write_bytes(b'changed executable')
        with self.assertRaisesRegex(ValueError, 'identical sealed solver'):
            nc.resume_inputs(self.job, 'midpoint', self.digest)

    def test_every_existing_numerical_and_physical_input_remains_bound(self):
        path = self.job/'inputs/case.json'; original = json.loads(path.read_text())
        for key, value in [('relative_force_tolerance', 1e-6), ('relative_linear_tolerance', 1e-8),
                           ('minimum_increment', 1e-9), ('tangent_perturbation', 1e-7)]:
            bad = copy.deepcopy(original); bad['settings'][key] = value
            path.write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError, 'input mismatch: stages'):
                nc.resume_inputs(self.job, 'midpoint', self.digest)
        for budget in (0., -1., 181., float('inf'), float('nan')):
            with self.assertRaises(AssertionError):
                nc.run(self.binary, self.folder/'invalid-budget', budget=budget)
        self.assertFalse((self.folder/'invalid-budget').exists())

    def test_omitted_optional_audit_and_changed_old_history(self):
        before = dict(attempts=[{'accepted': True, 'elapsed_seconds': 1.}], frames=[{'name': 'load'}],
                      geostatic={'attempts': []}, last_committed={'stage_progress': .5})
        after = copy.deepcopy(before); after['attempts'][0]['elapsed_seconds'] = 2.
        after['attempts'].append({'accepted': True}); after['last_committed']['stage_progress'] = .75
        row = nc.restart_observation(before, after)
        self.assertEqual(row['new_accepted'], 1)
        self.assertFalse(row['uninterrupted_full_path_compared'])
        after['attempts'][0]['accepted'] = False
        with self.assertRaisesRegex(ValueError, 'prefix changed: attempts'):
            nc.restart_observation(before, after)


if __name__ == '__main__':
    unittest.main()
