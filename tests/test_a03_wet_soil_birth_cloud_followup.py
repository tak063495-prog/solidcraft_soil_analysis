"""Small contracts and immutable evidence; no new FEM recomputation in CI."""
import hashlib
import importlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-12-30-a03-wet-soil-birth-cloud-followup'
sys.path.insert(0, str(HERE))
from a03follow_binding import bind, reference, fields
from a03follow_controls import fixture, controls, layout, dot, norm, common, compare
from a03follow_run import make_plan
from a03follow_report import effect_metric
from a03follow_surface import hydro_order


def read(name):
    return json.loads((HERE/name).read_text(encoding='utf-8'))


class A03CloudFollowup(unittest.TestCase):
    def test_formed_interaction_is_not_subtracted_twice(self):
        result = effect_metric(np.array([2.]), np.array([-4.]))
        self.assertEqual(result['absolute_l2'], 2.)
        self.assertEqual(result['relative_l2'], .5)

    def test_scope_and_null_accounting_are_retained(self):
        a = read('acceptance.json')
        for key in ('passed', 'a03_complete', 'new_limited_acceptance', 'production_rust_changed',
                    'nc_code_or_configuration_changed', 'acceptance_tolerance_changed', 'old_values_rewritten',
                    'physical_instability_proved', 'other_items_complete'):
            self.assertIs(a[key], False, key)
        self.assertIs(a['refusals_retained'], True)
        for value in read('measurements.json')['water_accounting'].values():
            self.assertIsNone(value)

    def test_document_terminal_status_with_crlf(self):
        for name in ('README.ja.md', 'REPRODUCE.ja.md'):
            s = (HERE/name).read_text(encoding='utf-8')
            self.assertTrue(s.endswith('A-03は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as d:
                p = Path(d)/name
                p.write_bytes(s.replace('\n', '\r\n').encode('utf-8'))
                self.assertEqual(p.read_text(encoding='utf-8'), s)

    def test_fixed_order_and_exact_common_coordinates(self):
        self.assertEqual(dot([1e16, 1., -1e16], [1., 1., 1.]), 1.)
        self.assertEqual(norm([[3., 0.], [0., 4.]]), 5.)
        rows = [(fields.row(s, i/8), {'u': np.array([i/8])}) for s in ('geostatic', 'fill', 'contact', 'reload') for i in range(1, 9)]
        selected = common(rows, 4)
        self.assertEqual(len(selected), 16)
        self.assertEqual(compare(selected, selected)['fields']['u']['maximum_absolute'], 0.)
        with self.assertRaises(ValueError):
            compare(selected, rows)

    def test_two_geometries_and_two_fresh_stress_distributions(self):
        sizes, stresses = set(), set()
        for shape in ('wide', 'tall'):
            for stress in ('k05', 'k08'):
                m, g, c = fixture(4, shape, stress)
                self.assertEqual((len(m['nodes']), len(m['tetrahedra10'])), (54, 12))
                sizes.add(tuple(np.max(m['nodes'], axis=0)))
                fresh = c['stages'][0]['states']['UpperClay']
                stresses.add(tuple(fresh['effective_stress_compression_positive_kpa']))
                self.assertEqual(fresh['plastic_strain_compression_positive'], [0.]*6)
                self.assertEqual(fresh['equivalent_plastic_strain'], 0.)
                self.assertEqual(fresh['last_return'], 'elastic')
                self.assertEqual(c['joint_system']['joints'][0]['water'], {'kind': 'continuous'})
                self.assertGreaterEqual(g['initialization']['water_level_m'], max(p[2] for p in m['nodes']))
        self.assertEqual((len(sizes), len(stresses)), (2, 2))

    def test_order_controls_restore_and_preserve_physical_keys(self):
        import fem_reference as fem
        saved = fem.volume, fem.face, fem.jacobi, reference.contact_work
        m, g, c = fixture(4, 'wide', 'k05')
        base = layout(reference.Model(m, g, c))
        with self.assertRaisesRegex(RuntimeError, 'retained'):
            with controls(7, True, 2):
                mm, gg, cc = fixture(4, 'wide', 'k05', True, True)
                model = reference.Model(mm, gg, cc)
                mapped = layout(model)
                self.assertEqual(base['element_keys'], mapped['element_keys'])
                self.assertEqual(base['contact_keys'], mapped['contact_keys'])
                self.assertEqual(len(model.ops[0][0]['weight']), 343)
                raise RuntimeError('retained')
        self.assertEqual(saved, (fem.volume, fem.face, fem.jacobi, reference.contact_work))
        import hydrostatic_reference as hydro
        rule = hydro.rule
        x, w = rule(4)
        with hydro_order(True):
            xx, ww = hydro.rule(4)
            np.testing.assert_array_equal(xx, x[::-1])
            np.testing.assert_array_equal(ww, w[::-1])
        self.assertIs(hydro.rule, rule)

    def test_full_factorial_and_more_common_points(self):
        plan = make_plan()['specs']
        self.assertEqual(len([s for s in plan if 'geometry' in s['roles']]), 12)
        self.assertEqual(len([s for s in plan if 'order' in s['roles']]), 8)
        cross = {(s['n'], s['minimum']) for s in plan if 'cross' in s['roles']}
        self.assertEqual(cross, {(n, m) for n in (4, 8, 16, 32) for m in (1, 2, 4)})
        data = read('measurements.json')
        self.assertEqual((data['formal_paths'], data['formal_states']), (60, 2688))
        self.assertEqual((data['formal_paths_including_full_surface'], data['formal_states_including_full_surface']), (64, 2816))
        self.assertTrue(all(p['process']['exit_code'] == 0 and not p['process']['timed_out'] for p in data['paths']))
        self.assertFalse(data['finest_grid_is_exact_reference'])

    def test_no_default_wiring_and_import_leaks(self):
        names = ('fields', 'reference', 'study')
        before = {n: sys.modules.get(n) for n in names}
        bind()
        self.assertEqual(before, {n: sys.modules.get(n) for n in names})
        for name in ('stages.rs', 'elements.rs', 'joint_law.rs'):
            self.assertNotIn('A03FOLLOW_', (ROOT/'src'/name).read_text(encoding='utf-8'))

    def test_existing_refusals_and_failed_values_are_unchanged(self):
        prior = read('prior-record-binding.json')
        for path, binding in prior['files'].items():
            self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(), binding['sha256'], path)

    def test_evidence_sha_without_equilibrium(self):
        index = read('archive-index.json')
        self.assertEqual(hashlib.sha256((HERE/'archive-index.json').read_bytes()).hexdigest(), read('acceptance.json')['archive_index_sha256'])
        for part in index['parts']:
            p = HERE/part['name']
            self.assertEqual(p.stat().st_size, part['size'])
            with p.open('rb') as f:
                self.assertEqual(hashlib.file_digest(f, 'sha256').hexdigest(), part['sha256'])
        self.assertIn('formal/difference-vectors.json', index['entries'])

    def test_source_access_is_explicit(self):
        review = read('primary-source-review.json')
        self.assertFalse(review['external_reference_values_introduced'])
        self.assertFalse(review['book_body_claims_invented'])
        self.assertTrue(any('Potts' in s['title'] and '本文未取得' in s['access'] for s in review['sources']))


if __name__ == '__main__':
    unittest.main()
