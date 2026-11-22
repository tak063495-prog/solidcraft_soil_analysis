"""Check a real nonlinear rejection and the subsequent reassembled equilibria."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parents[1] / 'interim/2026-11-22-impl-c01-c03-nonlinear'
spec = importlib.util.spec_from_file_location('nonlinear_c03_probe', HERE/'probe_c03.py')
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


class NonlinearJointProposals(unittest.TestCase):
    def test_failed_open_increment_reaches_same_anchor_with_unchanged_force_gate(self):
        old = study.measure(False, False)
        new = study.measure(False, True)
        self.assertIn('independent line search', old['stop'])
        self.assertEqual(old['rejected'][0], new['rejected'][0])
        self.assertEqual(old['prefix_sha256'], new['prefix_sha256'])
        self.assertEqual(new['stop'], 'observation point')
        self.assertEqual(new['last']['progress'], old['rejected'][0]['target'])
        self.assertGreater(len(new['observed_phase']), 1)
        self.assertLess(new['maximum_reassembled_force'], 1e-10)
        self.assertFalse(new['numerical_acceptance'])

    def test_mc_edge_rechecks_the_existing_principal_order_gate(self):
        import numpy as np
        import c03diag_mc as mc
        mat = study.read(HERE.parents[1] / 'examples/ground-joint-water-mc.json')['materials'][0]['mohr_coulomb']
        mat['cohesion_kpa'] = 30.
        # Strain from an admissible edge stress to a trial just outside both shear faces.
        bottom = (.5*82 - 2*30*np.cos(np.deg2rad(30))) / 1.5
        old = dict(effective_stress_compression_positive_kpa=[82., 82., bottom, 0., 0., 0.],
                   plastic_strain_compression_positive=[0.]*6, void_ratio=.7,
                   equivalent_plastic_strain=0., plastic_work_kpa=0.)
        nu = mat['poisson']
        compliance = (np.eye(3)*(1+nu) - nu) / mat['young_kpa']
        strain = np.r_[compliance @ [5e-6, 4.95e-6, -5e-6], [0.]*3]
        with self.assertRaisesRegex(ArithmeticError, 'inadmissible principal order'):
            mc.step(mat, old, strain)
        result = mc.step(mat, old, strain, require_order=True)
        stress = result['effective_stress_compression_positive_kpa']
        self.assertGreaterEqual(stress[0], stress[1] - 1e-8*(82+5e-6))
        self.assertLess(mc.normalized_yield(mat, stress), 1e-8)
        self.assertEqual(old['plastic_strain_compression_positive'], [0.]*6)

    def test_independent_proposals_resolve_all_anchors_and_check_both_half_steps(self):
        path = HERE.parent / '2026-12-12-c01-switch-event-localization/study.py'
        spec = importlib.util.spec_from_file_location('nonlinear_c01_study', path)
        c01 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(c01)
        examples = HERE.parents[1] / 'examples'
        mesh = c01.read(examples / 'joint-water-column-tet10-mesh.json')
        ground = c01.fem.normalize_ground(c01.read(examples / 'ground-joint-water-mc.json'))
        stages = c01.fem.normalize_stage_case(c01.read(examples / 'stages-joint-water-cycle.json'))
        steps = 8
        case = c01.prior.configuration(stages, steps, 5e-5)
        case['settings']['minimum_increment'] = min(case['settings']['minimum_increment'], 1/(steps*65536))
        control = c01.path_settings(True)
        case['joint_system']['global_step_control'] = dict(path_error=control)
        result = c01.reference._reference(mesh, ground, case, steps, 1, 7, False, True,
                                          maximum_m=5e-5, path_error=control)
        anchors = {(phase, i/steps) for phase in c01.STAGES for i in range(1, steps+1)}
        self.assertLessEqual(anchors, {(row['stage'], row['progress']) for row in result['path']})
        self.assertTrue(all(min(row['contact_partition_counts']) >= 4 for row in result['path']))
        self.assertTrue(any(not trial['accepted'] for trial in result['attempts']))
        for trial in result['attempts']:
            audit = trial['joint_path_error']
            self.assertGreater(audit['proposed_increment'], 0.)
            self.assertLess(max(trial['relative_force_residual'], audit['midpoint_relative_force_residual'],
                                audit['fine_relative_force_residual']), 1e-10)
            inside = (audit['maximum_slip_error_m'] <= audit['slip_budget_m']
                      and audit['maximum_dissipation_error_kpa_m'] <= audit['dissipation_budget_kpa_m'])
            self.assertEqual(trial['accepted'], inside)

    def test_documents_keep_unachieved_decision_under_windows_line_endings(self):
        for name in ('README.ja.md', 'REPRODUCE.ja.md'):
            original = (HERE/name).read_text(encoding='utf-8')
            for ending in ('\n', '\r\n'):
                with self.subTest(name=name, ending=repr(ending)), tempfile.TemporaryDirectory() as tmp:
                    path = Path(tmp)/name
                    path.write_bytes(original.replace('\n', ending).encode())
                    self.assertTrue(path.read_text(encoding='utf-8').endswith('C-01・C-03 は未達、passed=false。\n'))
