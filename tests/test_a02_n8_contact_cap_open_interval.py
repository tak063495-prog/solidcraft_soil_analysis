"""Small input/evidence contracts; the FEM paths are never rerun by CI."""
import importlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2027-01-05-a02-n8-contact-cap-open-interval'
sys.path.insert(0, str(HERE))
from a02n8_binding import norm, dot
from a02n8_run import schedule
from a02n8_evidence import verify


def read(name):
    return json.loads((HERE/name).read_text(encoding='utf-8'))


class A02N8ContactCapOpenInterval(unittest.TestCase):
    def test_incomplete_scope_and_preserved_guards(self):
        a = read('acceptance.json')
        for key in ('passed', 'a02_complete', 'new_limited_acceptance', 'production_rust_changed',
                    'nc_code_or_configuration_changed', 'acceptance_tolerance_changed', 'old_values_rewritten',
                    'other_items_complete', 'exact_continuous_switch_roots_identified'):
            self.assertIs(a[key], False, key)
        self.assertIs(a['refusals_retained'], True)

    def test_prescribed_grid_only_subdivides_first_open_interval(self):
        for n in (4, 8):
            original = [i/n for i in range(1, n+1)]
            for b in (2, 4, 8):
                grid = schedule(n, b)
                self.assertEqual(set(grid), {'open'})
                self.assertEqual(grid['open'][b-1:], original)
                self.assertEqual(grid['open'][:b], [i/(n*b) for i in range(1, b+1)])
                self.assertTrue(all(a < c for a, c in zip(grid['open'], grid['open'][1:])))
        with self.assertRaises(ValueError):
            schedule(768, 8)

    def test_fixed_order_observations(self):
        self.assertEqual(norm([[3., 0.], [0., 4.]]), 5.)
        self.assertEqual(dot([1e16, 1., -1e16], [1., 1., 1.]), 1.)
        self.assertEqual(dot([[1., 2.], [3., 4.]], [[4., 3.], [2., 1.]]), 20.)

    def test_all_requested_exchanges_are_sealed_observations(self):
        m = read('measurements.json')
        paths = m['paths']
        exchange = [p for p in paths if p['path'].startswith('exchange/')]
        self.assertEqual(len(exchange), 16)
        self.assertEqual(len({(p['cap'], p['engine'], p['predictor'], p['repartition']) for p in exchange}), 16)
        self.assertTrue(all(p['n'] == 8 and p['parts'] == 1 and p['complete'] and p['exit_code'] == 0 for p in exchange))
        self.assertEqual(len(m['exchange']), 8)
        for row in m['exchange'].values():
            points = row['cap_trace']['accepted_points']
            self.assertEqual(len(points), 56)
            self.assertTrue(all(len(p['contact_parts']) == 32 for p in points))
            self.assertGreater(row['cap_trace']['maximum_contact_parts'], 1)
            self.assertEqual(row['switches_cap_vs_normal']['only_left'], 0)
            self.assertEqual(row['switches_cap_vs_normal']['only_right'], 0)
        deferred = m['exchange']['native-pnative-rindependent']['cap_trace']
        self.assertGreater(deferred['repartitions'], 0)
        self.assertEqual(sum(p['emitted_repartition_events'] for p in deferred['accepted_points']), 0)

    def test_open_bisections_are_fresh_equilibrium_paths(self):
        m = read('measurements.json')
        paths = [p for p in m['paths'] if p['path'].startswith('open/')]
        self.assertEqual({(p['n'], p['divisions'], p['engine']) for p in paths},
                         {(n, b, e) for n in (4, 8) for b in (2, 4, 8) for e in ('native', 'independent')})
        for p in paths:
            self.assertTrue(p['complete'])
            self.assertEqual(p['accepted_states'], 7*p['n']+p['divisions']-1)
            self.assertEqual(p['parts'], 4)
        self.assertEqual(m['principal_accepted_states'], 1444)
        self.assertEqual(m['observer_free_accepted_states'], 722)
        self.assertIs(m['observations_are_bracket_upper_endpoints_not_exact_continuous_roots'], True)
        for engine in ('native', 'independent'):
            finest = m['open']['open/n8-b8/'+engine]['observations']
            self.assertEqual(finest['branch_summary']['unique_slip_activity_gp'], 32)
            self.assertEqual(finest['branch_summary']['slip_activity_interval_gp_occurrences'], 39)
            active = [r for r in finest['intervals'] if r['slip_activity_gp']]
            self.assertEqual([(r['lower'], r['upper'], r['slip_activity_gp'], r['opening_gp']) for r in active],
                             [(.03125, .046875, 32, 25), (.046875, .0625, 7, 7)])

    def test_document_tails_with_real_crlf(self):
        for name in ('README.ja.md', 'REPRODUCE.ja.md'):
            text = (HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('A-02は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as folder:
                path = Path(folder)/name
                path.write_bytes(text.replace('\n', '\r\n').encode('utf-8'))
                self.assertTrue(path.read_text(encoding='utf-8').endswith('A-02は未達、passed=false。\n'))

    def test_no_production_hook_or_generic_import_leak(self):
        for path in (ROOT/'src').rglob('*.rs'):
            self.assertNotIn('A02N8_', path.read_text(encoding='utf-8'))
        names = ('reference', 'fields', 'study', 'projection', 'geometry_audit')
        before = {k: sys.modules.get(k) for k in names}
        importlib.reload(importlib.import_module('a02n8_binding'))
        self.assertEqual(before, {k: sys.modules.get(k) for k in names})

    def test_archive_and_primary_source_limits(self):
        self.assertIs(verify(HERE)['all_sha_verified'], True)
        review = read('primary-source-review.json')
        self.assertIs(review['external_reference_values_introduced'], False)
        for title in ('Wriggers', 'Laursen'):
            self.assertTrue(any(title in s['title'] and '本文未取得' in s['access'] for s in review['sources']))


if __name__ == '__main__':
    unittest.main()
