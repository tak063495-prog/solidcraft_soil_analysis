"""Explicit opt-in tests; outside existing tests/tools and pytest test_*.py names."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from reference import (CASES, contact_return, elasticity, face_rule, make_model,
                       run, simplex, dump)
from compare import (Refusal, compare, export_csv, field_differences, load,
                     validate, remap, TOLERANCES)
from study import fields, differences, summarize, step_difference, gp_fields
from audit import audit


class Algebra(unittest.TestCase):
    def test_tri3_affine_strain_and_energy(self):
        x = np.array([[0., 0.], [2., 0.], [0., 1.]])
        b, volume = simplex(x)
        u = np.column_stack((.125*x[:, 0]+.25*x[:, 1], -.125*x[:, 1]))
        strain = b@u.ravel()
        self.assertTrue(np.array_equal(strain, [.125, -.125, .25]))
        self.assertEqual(volume, 1.)
        self.assertEqual(.5*volume*strain@elasticity(2)@strain, 75.)

    def test_tet4_affine_strain_and_volume(self):
        x = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        b, volume = simplex(x)
        strain = b@(.125*x).ravel()
        self.assertTrue(np.array_equal(strain, [.125, .125, .125, 0., 0., 0.]))
        self.assertEqual(volume, 1.)
        self.assertEqual(.5*strain@elasticity(3)@strain, 140.625)

    def test_rigid_translation_strain_zero(self):
        for dim in (2, 3):
            x = np.vstack((np.zeros(dim), np.eye(dim)))
            b, _ = simplex(x)
            self.assertTrue(np.array_equal(b@np.ones(dim*(dim+1)), np.zeros(len(b))))

    def test_segment_quadrature_exact_measure(self):
        for order in (1, 2, 3):
            rule = face_rule(np.array([[0., 0.], [2., 0.]]), order)
            self.assertLessEqual(abs(sum(w for _, w in rule)-2), TOLERANCES['independent_elastic_relative'])

    def test_tri3_quadrature_integrates_shape_at_order_two(self):
        rule = face_rule(np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]), 2)
        value = sum(n*w for n, w in rule)
        self.assertLessEqual(np.linalg.norm(value-np.ones(3)/6), TOLERANCES['independent_elastic_relative'])

    def test_closed_form_stick_slip_reversal(self):
        n = np.array([1., 0.]); old = np.zeros(2)
        q, _, slip, state, energy = contact_return(np.array([-.125, .0625]), old, n, 8., 4., .5)
        self.assertEqual(state, 'stick'); self.assertTrue(np.array_equal(q, [-1., .25]))
        q, _, slip, state, _ = contact_return(np.array([-.125, .25]), old, n, 8., 4., .5)
        self.assertEqual(state, 'slip'); self.assertTrue(np.array_equal(q, [-1., .5]))
        self.assertTrue(np.array_equal(slip, [0., .125]))
        q, _, slip, state, _ = contact_return(np.array([-.125, -.25]), slip, n, 8., 4., .5)
        self.assertTrue(np.array_equal(q, [-1., -.5])); self.assertEqual(state, 'slip')
        self.assertTrue(np.array_equal(slip, [0., -.125]))

    def test_open_contact_is_unilateral_and_does_not_commit_trial_history(self):
        old = np.array([0., .125]); saved = old.copy()
        q, k, slip, state, energy = contact_return(np.array([.125, 1.]), old, np.array([1., 0.]), 8., 4., .5)
        self.assertTrue(np.array_equal(q, np.zeros(2))); self.assertTrue(np.array_equal(k, np.zeros((2, 2))))
        self.assertEqual(energy, 0.); self.assertEqual(state, 'open')
        self.assertTrue(np.array_equal(old, saved)); self.assertTrue(np.array_equal(slip, old))

    def test_consistent_coulomb_tangent_by_directional_derivative(self):
        n = np.array([1., 0., 0.]); old = np.zeros(3)
        jump = np.array([-.125, .25, .125]); delta = np.array([.25, -.125, .25])
        q, k, *_ = contact_return(jump, old, n, 8., 4., .5)
        h = 2.**-20
        qp = contact_return(jump+h*delta, old, n, 8., 4., .5)[0]
        qm = contact_return(jump-h*delta, old, n, 8., 4., .5)[0]
        self.assertLess(np.linalg.norm((qp-qm)/(2*h)-k@delta), TOLERANCES['independent_elastic_relative'])

    def test_shared_junction_has_one_physical_node(self):
        model = make_model('shared2d')
        self.assertEqual(sum(np.array_equal(x, [1., 1.]) for x in model.x), 1)
        self.assertEqual(len(model.node_ids), 33)
        u = np.ones_like(model.x)
        for p in model.contacts:
            self.assertTrue(np.array_equal(p['shape']@(u[p['plus']]-u[p['minus']]), [0., 0.]))

    def test_axial_member_formula_and_end_force_moment(self):
        model = make_model('beam2d')
        m = model.members[0]
        displacement = np.array([0., 0., .125, 0.])
        force = m['ea']/m['length']*(m['b']@displacement)
        self.assertEqual(force, 75.)
        self.assertTrue(np.array_equal(m['b']*force, [-75., 0., 75., 0.]))


class Histories(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.docs = {name: run(name, 4) for name in CASES}

    def test_all_cases_complete_with_slip_and_no_expanded_acceptance(self):
        for name, doc in self.docs.items():
            with self.subTest(case=name):
                self.assertTrue(doc['complete']); self.assertIsNone(doc['accepted'])
                self.assertGreaterEqual(len(doc['joint_ids']), 2)
                self.assertGreater(max(r['accumulated_slip_m'] for r in doc['frames'][-1]['contacts']), 0.)
                self.assertEqual(len(doc['frames']), 12)
                validate(doc)

    def test_every_nodal_force_and_independent_stress_scatter(self):
        for doc in self.docs.values():
            report = summarize(doc)
            self.assertLessEqual(report['maximum_node_residual_kn'], TOLERANCES['nodal_force_kn'])
            self.assertLessEqual(report['maximum_matrix_vs_stress_scatter_kn'], TOLERANCES['nodal_force_kn'])

    def test_tensor_coulomb_axial_and_discrete_work_audit(self):
        for doc in self.docs.values():
            report = audit(doc)
            self.assertTrue(report['complete']); self.assertIsNone(report['accepted'])
            self.assertFalse(report['physical_work_corrected'])
            self.assertEqual(report['state_mismatches'], 0)
            for key in ('bulk_tensor_scatter_kn', 'contact_scatter_kn', 'section_force_kn'):
                self.assertLessEqual(report['maxima'][key], TOLERANCES['nodal_force_kn'])
            for key in ('friction_work_identity_kn_m', 'closed_discrete_work_identity_kn_m',
                        'member_work_identity_kn_m', 'global_work_accounting_kn_m'):
                self.assertLessEqual(report['maxima'][key], TOLERANCES['work_kn_m'])

    def test_elastic_preload_work_identity_uses_quoted_work_diagnostic(self):
        for doc in self.docs.values():
            for frame in doc['frames'][:4]:
                self.assertLessEqual(abs(frame['diagnostics']['work_defect_kn_m']), TOLERANCES['work_kn_m'])
                self.assertEqual(frame['friction_dissipation_kn_m'], 0.)

    def test_sliding_work_defect_is_retained_not_implicitly_corrected(self):
        for doc in self.docs.values():
            f = doc['frames'][-1]
            energy = sum(f[k] for k in ('bulk_energy_kn_m', 'contact_energy_kn_m', 'member_energy_kn_m'))
            self.assertEqual(f['diagnostics']['work_defect_kn_m'], f['external_work_kn_m']-energy-f['friction_dissipation_kn_m'])
            self.assertFalse(summarize(doc)['quoted_work_diagnostic'])
            self.assertEqual(doc['transformations'], [])

    def test_beam_and_anchor_share_load_and_change_contact_solution(self):
        for name, baseline in [('beam2d', 'parallel2d'), ('anchor3d', 'parallel3d')]:
            f = self.docs[name]['frames'][-1]
            self.assertTrue(all(m['axial_force_kn'] != 0. for m in f['members']))
            self.assertNotEqual(f['nodes'], self.docs[baseline]['frames'][-1]['nodes'])
            self.assertGreater(f['member_energy_kn_m'], 0.)

    def test_contact_joint_node_order_raw_deltas_are_recorded(self):
        base = self.docs['cross2d']; a, states = fields(base)
        for variant in ('contacts', 'joints', 'nodes'):
            doc = run('cross2d', 4, variant=variant); b, changed_states = fields(doc)
            metrics = differences(a, b)
            # Legacy relative diagnostic is reported; no new slip/order criterion.
            self.assertLess(np.linalg.norm(a['displacement_m']-b['displacement_m'])/np.linalg.norm(a['displacement_m']),
                            TOLERANCES['independent_elastic_relative'])
            self.assertEqual(states, changed_states)
            self.assertTrue(all(v['accepted'] is None and v['tolerance'] is None for v in metrics.values()))

    def test_step_and_gp_series_never_claim_pointwise_or_numerical_acceptance(self):
        fine = run('parallel2d', 8)
        result = step_difference(self.docs['parallel2d'], fine)
        self.assertEqual(result['common_actual_endpoints'], 12); self.assertIsNone(result['accepted'])
        coarse_gp = run('parallel2d', 4, 1)
        delta = differences(gp_fields(coarse_gp), gp_fields(self.docs['parallel2d']))
        self.assertNotIn('traction_kpa', delta)
        self.assertIsNone(delta['joint_resultant_kn']['tolerance'])


class ExportContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = run('beam2d', 4)
        cls.mapping = {jid: jid for jid in cls.doc['joint_ids']}

    def witness(self, doc=None):
        value = copy.deepcopy(self.doc if doc is None else doc)
        value['evaluation_order'].reverse()
        for f in value['frames']:
            f['contacts'].reverse(); f['nodes'].reverse(); f['members'].reverse()
        return value

    def test_csv_and_json_roundtrip_preserve_all_raw_fields(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'raw.csv'; export_csv(self.doc, p)
            actual = load(p, p.with_suffix('.manifest.json'))
            self.assertEqual(actual, self.doc)
            dump(Path(folder)/'raw.json', actual)
            self.assertEqual(load(Path(folder)/'raw.json'), self.doc)

    def test_reordered_serialization_is_identity_but_not_physical_acceptance(self):
        result = compare(self.doc, self.doc, self.mapping, self.witness())
        self.assertTrue(result['exact_raw_match']); self.assertTrue(result['exact_order_match'])
        self.assertFalse(result['comparison_refused']); self.assertIsNone(result['accepted'])

    def test_units_signs_corrections_and_incomplete_output_are_refused(self):
        mutations = [lambda d: d['units'].update(force='N'),
                     lambda d: d['signs'].update(jump='u_minus-u_plus'),
                     lambda d: d['transformations'].append('clip'),
                     lambda d: d.update(raw=False), lambda d: d.update(complete=False),
                     lambda d: d.pop('producer'), lambda d: d['evaluation_order'].pop()]
        for mutate in mutations:
            bad = copy.deepcopy(self.doc); mutate(bad)
            with self.assertRaises(Refusal): validate(bad)

    def test_missing_duplicate_joint_gp_generation_node_member_fields_are_refused(self):
        mutations = [lambda d: d['frames'][0]['contacts'].pop(),
                     lambda d: d['frames'][0]['contacts'].append(copy.deepcopy(d['frames'][0]['contacts'][0])),
                     lambda d: d['frames'][0]['contacts'][0].pop('generation'),
                     lambda d: d['frames'][0]['nodes'].pop(),
                     lambda d: d['frames'][0]['members'].pop(),
                     lambda d: d['frames'][0]['contacts'][0].update(joint_id='unknown'),
                     lambda d: d['frames'][0]['contacts'][0].update(traction_kpa=[float('nan'), 0.])]
        for mutate in mutations:
            bad = copy.deepcopy(self.doc); mutate(bad)
            with self.assertRaises((Refusal, KeyError)): validate(bad)

    def test_joint_id_mapping_is_explicit_bijective_and_complete(self):
        changed = remap(self.doc, {jid: 'export-'+jid for jid in self.doc['joint_ids']})
        mapping = {'export-'+jid: jid for jid in self.doc['joint_ids']}
        result = compare(self.doc, changed, mapping, self.witness(changed))
        self.assertTrue(result['exact_raw_match'])
        with self.assertRaises(Refusal): compare(self.doc, changed, {}, self.witness(changed))
        with self.assertRaises(Refusal): remap(changed, {k: 'same' for k in mapping})

    def test_numerical_smoothing_even_with_false_raw_declaration_is_not_accepted(self):
        changed = copy.deepcopy(self.doc)
        changed['frames'][-1]['contacts'][0]['plastic_slip_m'][1] += .001
        result = compare(self.doc, changed, self.mapping, self.witness(changed))
        self.assertTrue(result['comparison_refused']); self.assertTrue(result['differences'])
        self.assertIsNone(result['accepted'])

    def test_order_dependent_result_is_refused_and_witness_is_required(self):
        with self.assertRaises(Refusal): compare(self.doc, self.doc, self.mapping, None)
        with self.assertRaises(Refusal): compare(self.doc, self.doc, self.mapping, self.doc)
        changed = self.witness(); changed['frames'][-1]['members'][0]['axial_force_kn'] += 1.
        result = compare(self.doc, self.doc, self.mapping, changed)
        self.assertTrue(result['comparison_refused']); self.assertTrue(result['order_differences'])

    def test_csv_corruption_hash_is_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'raw.csv'; export_csv(self.doc, p)
            p.write_bytes(p.read_bytes()+b'\n')
            with self.assertRaises(Refusal): load(p, p.with_suffix('.manifest.json'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
