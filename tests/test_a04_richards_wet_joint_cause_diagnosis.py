"""Small observation contracts; full causal solves remain sealed outside CI."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis'
sys.path.insert(0, str(HERE))
from a04cause_binding import bind
from a04cause_metrics import common, compare, convergence, metric, norm
from a04cause_driver import reverse_keys
from a04cause_report import canonical, control_rows


class A04RichardsWetJointCauseDiagnosis(unittest.TestCase):
    def test_null_or_nonfinite_is_not_a_zero_water_budget(self):
        for v in (None, math.nan, math.inf):
            with self.assertRaises(ValueError):
                metric([v], [0.])
        with self.assertRaises(ValueError):
            metric([1.], [[1.]])

    def test_fixed_order_norm_and_original_denominator(self):
        self.assertEqual(norm([[3., 4.], [0., 12.]]), 13.)
        v = metric([1e-18], [0.])
        self.assertEqual(v['reference_l2'], 0.)
        self.assertEqual(v['denominator'], 1e-12)
        self.assertEqual(v['maximum_absolute'], 1e-18)
        self.assertEqual(v['relative_l2'], 1e-18/1e-12)
        self.assertTrue(v['denominator_floor_active'])

    def test_physical_time_and_release_progress_are_distinct(self):
        row = lambda t, p: (dict(phase='release', time_s=t, progress=p), dict(u=[p]))
        base = [row(2., .5), row(2., 1.)]
        other = [row(2., .25), *base]
        self.assertEqual(compare(base, common(other, base))['states'], 2)
        with self.assertRaises(ValueError):
            common([row(.5, .5), row(1., 1.)], base)
        with self.assertRaises(ValueError):
            common([*base, base[0]], base)

    def test_no_interpolation_or_missing_ID_field(self):
        a = [(dict(time_s=.25), dict(cell_ids=[1], u=[[0., 1.]]))]
        b = [(dict(time_s=.125), a[0][1])]
        with self.assertRaises(ValueError):
            compare(a, b)
        with self.assertRaises(ValueError):
            compare(a, [(a[0][0], dict(u=[[0., 1.]]))])
        c = [(a[0][0], dict(cell_ids=[2], u=[[0., 1.]]))]
        self.assertEqual(compare(a, c)['fields']['cell_ids']['mismatches'], 1)

    def test_quadrature_control_compares_only_shared_unknowns(self):
        row = dict(time_s=1.)
        a = [(row, dict(u=[1.], hydraulic_positions=[[0., 0., 0.]]))]
        b = [(row, dict(u=[2.], hydraulic_positions=[[0., 0., 0.], [1., 1., 1.]]))]
        self.assertEqual(compare(a, b, ['u'])['fields']['u']['maximum_absolute'], 1.)
        with self.assertRaises(ValueError):
            compare(a, b)

    def test_finite_series_and_zero_difference_do_not_make_an_exact_reference(self):
        grid = [(dict(time_s=.5), dict(u=[1.]))]
        data = convergence([grid, grid, grid], grid, ['u'])
        self.assertIsNone(data['orders']['u'])
        self.assertFalse(data['numerical_acceptance'])
        self.assertFalse(data['extrapolation'])

    def test_generic_imports_and_search_path_are_restored(self):
        names = ('soil_reference', 'study', 'fields', 'reference')
        before = {n:sys.modules.get(n) for n in names}
        path = sys.path[:]
        modules = bind()
        self.assertEqual(set(modules), {'soil_reference', 'study'})
        for n in names:
            self.assertIs(sys.modules.get(n), before[n])
        self.assertEqual(sys.path, path)

    def test_reproduction_comparison_excludes_only_elapsed(self):
        a = dict(time_s=2., generation=1, head=-0., elapsed_seconds=1., state=dict(active=True))
        self.assertEqual(canonical(a), canonical(reverse_keys({**a, 'elapsed_seconds':9.})))
        self.assertNotEqual(canonical(a), canonical({**a, 'generation':2}))
        self.assertNotEqual(canonical(a), canonical({**a, 'head':0.}))
        self.assertNotEqual(canonical(a), canonical({**a, 'state':dict(active=1)}))

    def test_control_rows_keep_progress_but_do_not_create_physical_seconds(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'rows.jsonl'
            source = dict(row=dict(stage='geostatic',progress=.25,stage_index=0,factor=-.001,
                                   relative_force_residual=1e-17,contact_partitions=[1]), fields=dict(u=[0.]))
            p.write_text(json.dumps(source)+'\n',encoding='utf-8')
            self.assertEqual(control_rows(p)[0][0],dict(stage='geostatic',progress=.25))
            source['row']['time_s'] = .25
            p.write_text(json.dumps(source)+'\n',encoding='utf-8')
            with self.assertRaises(ValueError):
                control_rows(p)

    def test_document_terminal_status_with_real_CRLF(self):
        text = (HERE/'README.ja.md').read_text(encoding='utf-8')
        self.assertTrue(text.endswith('A-04は未達、passed=false。\n'))
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'README.ja.md'
            p.write_bytes(text.replace('\n', '\r\n').encode('utf-8'))
            self.assertEqual(p.read_text(encoding='utf-8'), text)

    def test_sealed_evidence_and_unachieved_interface_contracts(self):
        a = json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
        self.assertFalse(a['passed'])
        self.assertFalse(a['a04_complete'])
        index = HERE/'archive-index.json'
        self.assertEqual(hashlib.sha256(index.read_bytes()).hexdigest(), a['archive_index_sha256'])
        m = json.loads((HERE/'measurements.json').read_text(encoding='utf-8'))
        self.assertFalse(m['all_field_numerical_acceptance'])
        for k, v in m['interface'].items():
            if k == 'bidirectional_controller':
                self.assertFalse(v)
            else:
                self.assertIsNone(v)
        self.assertEqual(len(m['scope_refusals']['actual_cli_refusals']), 16)
        self.assertTrue(m['same_surface_boundary_gap']['actual_production_complete'])
        self.assertFalse(m['same_surface_boundary_gap']['actual_production_refused'])
        data = json.loads(index.read_text(encoding='utf-8'))
        for part in data['parts']:
            p = HERE/part['name']
            self.assertEqual(p.stat().st_size, part['size'])
            with p.open('rb') as f:
                self.assertEqual(hashlib.file_digest(f,'sha256').hexdigest(), part['sha256'])
        for required in ('formal/soil/VG-n16/native-fields.jsonl',
                         'formal/VG-quadrature7/independent-fields.jsonl',
                         'formal/boundary-gap/native.json'):
            self.assertIn(required, data['entries'])
