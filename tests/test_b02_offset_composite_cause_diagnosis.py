"""Small scientific/portability contracts; sealed solves are not repeated in CI."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-12-24-b02-offset-composite-cause-diagnosis'
sys.path.insert(0, str(HERE))
from b02cause_binding import bind, legacy
from b02cause_driver import reverse_keys
from b02cause_metrics import common, document_text, dot, metric, norm
from b02cause_reference import volume_control
from b02cause_report import holding_vector, probe_metric, structure_norms


class B02OffsetCompositeCauseDiagnosis(unittest.TestCase):
    def test_fixed_order_norm_retains_denominator_and_rejects_invalid_data(self):
        self.assertEqual(dot([1e16, 1., -1e16], [1., 1., 1.]), 1.)
        self.assertEqual(norm([[3., 4.], [0., 12.]]), 13.)
        value = metric([1e-18], [0.])
        self.assertEqual(value['relative_l2'], 1e-18/1e-12)
        self.assertEqual(value['reference_l2'], 0.)
        self.assertTrue(value['denominator_floor_active'])
        with self.assertRaises(ValueError):
            metric([math.nan], [0.])
        with self.assertRaises(ValueError):
            metric([0.], [0., 1.])

    def test_common_coordinates_are_exact_and_unique(self):
        rows = [dict(stage_index=s, progress=j/8) for s in range(5) for j in range(1, 9)]
        self.assertEqual(len(common(rows)), 20)
        with self.assertRaises(ValueError):
            common(rows+[rows[-1]])
        with self.assertRaises(ValueError):
            common(rows[:-1])
        probes = {(0, 1.): {'q': []}, (1, 1.): {'q': [3., 4.]}}
        self.assertEqual(probe_metric(probes, probes)['absolute_l2'], 0.)
        self.assertEqual(probe_metric(probes, probes)['components'], 2)
        with self.assertRaises(ValueError):
            probe_metric(probes, {(0, 1.): {'q': [0.]}, (1, 1.): {'q': [3., 4.]}})

    def test_import_names_and_sys_path_are_restored(self):
        names = ('cases', 'compare', 'reference', 'study', 'work_audit')
        before = {n: sys.modules.get(n) for n in names}
        path = sys.path[:]
        self.assertEqual(set(bind()), set(names))
        self.assertEqual(path, sys.path)
        for name in names:
            self.assertIs(before[name], sys.modules.get(name))
        obj = {'a': {'x': 1, 'y': 2}, 'b': [1, 2]}
        self.assertEqual(obj, reverse_keys(obj))
        self.assertEqual(list(reverse_keys(obj)), ['b', 'a'])

    def test_quadrature_control_is_local_and_restored_on_error(self):
        reference = legacy['reference']
        original = reference.ref
        mesh, _, _, _ = legacy['cases'].inputs()
        with self.assertRaisesRegex(RuntimeError, 'retained'):
            with volume_control(7, True):
                q = reference.ref.volume(mesh, 0, independent=True)
                self.assertEqual(len(q['weight']), 343)
                self.assertIs(legacy['cases'].ref, original)
                raise RuntimeError('retained')
        self.assertIs(reference.ref, original)

    def test_signed_offset_coupling_and_load_work_conjugacy(self):
        self.assertEqual(structure_norms([3., 4., 0., 1., 2., 2.]), {'force_kn': 5., 'moment_kn_m': 3.})
        np.testing.assert_array_equal(holding_vector({}, 3), [0., 0., 0.])
        np.testing.assert_array_equal(holding_vector({'transition_release_load_kn': []}, 3), [0., 0., 0.])
        np.testing.assert_array_equal(holding_vector({'transition_release_load_kn': [1., -2., 3.]}, 3), [1., -2., 3.])
        with self.assertRaises(ValueError):
            holding_vector({'transition_release_load_kn': [1.]}, 3)
        mesh, _, case, _ = legacy['cases'].inputs()
        kernel = legacy['reference'].System(mesh, case['beam_system']).plate_data[0][2]
        d = np.arange(36, dtype=float)*1e-7
        neutral = kernel.map@d
        force = kernel.neutral.load([.2, -.1, -.1])
        self.assertAlmostEqual(dot(kernel.map.T@force, d), dot(force, neutral), places=18)
        local_s = kernel.neutral.transform@d
        local_p = kernel.neutral.transform@neutral
        bm, bb, _, _ = kernel.neutral.operators([.2, .3, .5])
        np.testing.assert_allclose(bm@local_p, bm@local_s-kernel.e*(bb@local_p), atol=1e-20, rtol=1e-12)

    def test_document_suffix_contract_accepts_actual_crlf(self):
        marker = 'B-02は未達、passed=false。\n'
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'README.ja.md'
            p.write_bytes(('記録\r\n'+marker.replace('\n', '\r\n')).encode())
            self.assertTrue(document_text(p).endswith(marker))
        self.assertTrue(document_text(HERE/'README.ja.md').endswith(marker))

    def test_acceptance_and_original_failures_stay_distinct(self):
        value = json.loads((HERE/'acceptance.json').read_text())
        self.assertFalse(value['passed'])
        self.assertFalse(value['b02_complete'])
        self.assertIsNone(value['unsupported_combined_equilibrium'])
        self.assertEqual(value['new_maximum_tet10'], 96)
        old = ROOT/'interim/2026-12-01-tri6-offset-composite-b02/measurements.json'
        binding = json.loads((HERE/'prior-record-binding.json').read_text())
        self.assertEqual(hashlib.sha256(old.read_bytes()).hexdigest(), binding[old.relative_to(ROOT).as_posix()])
        for p in HERE.rglob('*.json'):
            if p.name in ('acceptance.json', 'measurements.json', 'validation.json'):
                self.assertEqual(p.parent, HERE, 'retained copies must stay sealed or be uniquely named for P6 inventory')

    def test_sealed_evidence_parts_exist_and_match_sha256(self):
        value = json.loads((HERE/'archive-index.json').read_text())
        h = hashlib.sha256()
        size = 0
        for part in value['parts']:
            raw = (HERE/part['name']).read_bytes()
            self.assertEqual(len(raw), part['size'])
            self.assertEqual(hashlib.sha256(raw).hexdigest(), part['sha256'])
            h.update(raw)
            size += len(raw)
        self.assertEqual(h.hexdigest(), value['archive_sha256'])
        self.assertEqual(size, value['archive_bytes'])


if __name__ == '__main__':
    unittest.main()
