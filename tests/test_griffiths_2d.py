"""Real archived displacement regression plus synthetic comparator controls.

No Griffiths FEM solve is invoked, including by the series preflight tests.
"""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import griffiths_2d as comparison
import prepare_slope_h_series as preparation
from plastic_spatial_convergence import read, write, sha

REFERENCE = ROOT/'tests/data/griffiths_2d_reference.json'


class PinnedOriginal(unittest.TestCase):
    def test_original_661_nodes_and_all_retained_copies_are_identical(self):
        ref = comparison.load_reference(REFERENCE)
        self.assertEqual([row[0] for row in ref['nodes']], list(range(1, 662)))
        self.assertEqual(len(ref['displacement_xz_m']), 661)
        self.assertEqual(ref['material_reduction_factor'], 1.3500000000000003)
        self.assertEqual(ref['source']['result_sha256'],
                         '8d603e36d17aff6fd21dd46adefe40a778364bcbe6c8eab1608860ea835e70ed')
        target = np.asarray(ref['displacement_xz_m'])
        self.assertTrue(np.isfinite(target).all())
        self.assertGreater(np.linalg.norm(target), 0)
        self.assertTrue(np.all(target[np.asarray(ref['nodes'])[:, 3:5].astype(bool)] == 0))
        copies = [ROOT/name for name in comparison.ARCHIVES]
        copies.append(ROOT/'interim/2026-10-30-Root635-full-audits-slope-reference/slope/coarse-field-samples.json')
        for path in copies:
            self.assertEqual(read(path)['source2d_u_xz_m'], ref['displacement_xz_m'])

    def test_reference_missing_or_any_content_changed_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'reference file required'):
            comparison.load_reference(None)
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'ref.json'
            with self.assertRaises(FileNotFoundError):
                comparison.load_reference(path)
            for change in ('value', 'count', 'node_order', 'units', 'factor', 'source'):
                ref = read(REFERENCE)
                if change == 'value': ref['displacement_xz_m'][42][0] += 1e-15
                if change == 'count': ref['displacement_xz_m'].pop()
                if change == 'node_order': ref['nodes'].reverse()
                if change == 'units': ref['units'] = 'mm'
                if change == 'factor': ref['material_reduction_factor'] = 1.
                if change == 'source': ref['source']['result_sha256'] = 'unverified'
                path.write_text(json.dumps(ref), encoding='utf-8')
                with self.subTest(change=change), self.assertRaisesRegex(ValueError, 'changed or unverified'):
                    comparison.load_reference(path)

    def test_reference_crlf_checkout_is_accepted(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'ref.json'
            path.write_bytes(REFERENCE.read_bytes().replace(b'\r\n', b'\n').replace(b'\n', b'\r\n'))
            self.assertEqual(comparison.load_reference(path), read(REFERENCE))

    def test_actual_archive_measurements_preserve_four_failures_and_indeterminate_continuum(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'actual.json'
            report = comparison.compare_archives(REFERENCE, output)
            self.assertTrue(report['comparison_measured'])
            self.assertFalse(report['new_slope_solve_performed'])
            self.assertFalse(report['griffiths_slope_five_percent_closed'])
            self.assertFalse(report['p6_complete'])
            self.assertIsNone(report['fs']['value'])
            self.assertIsNone(report['fs']['failure_upper_bound'])
            self.assertEqual(report['continuum_verdict'], 'indeterminate')
            for case, percent in zip(report['cases'], [15.58199870691725, 15.885691929516224,
                                                       15.593219845164077, 15.672391195315055], strict=True):
                self.assertEqual(case['five_percent_verdict'], 'fail')
                self.assertFalse(case['in_plane_h_series'])
                self.assertAlmostEqual(case['displacement_sample_relative_L2_per_depth']['0.5']*100, percent, places=11)
            with self.assertRaises(FileExistsError):
                comparison.compare_archives(REFERENCE, output)

    def test_archived_native_field_modification_is_refused(self):
        original = comparison.read
        def modified(path):
            data = original(path)
            if 'native3d_u_xz_m' in data:
                data['native3d_u_xz_m']['0.5'][42][0] += .001
            return data
        with tempfile.TemporaryDirectory() as temp, patch.object(comparison, 'read', side_effect=modified):
            out = Path(temp)/'no.json'
            with self.assertRaisesRegex(ValueError, 'archived displacement data changed'):
                comparison.compare_archives(REFERENCE, out)
            self.assertFalse(out.exists())

    def test_cli_requires_explicit_data_and_returns_failure_for_actual_archives(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'report.json'
            cmd = [sys.executable, str(ROOT/'tools/griffiths_2d.py'), 'archives', str(output)]
            self.assertEqual(subprocess.run(cmd, capture_output=True).returncode, 2)
            self.assertFalse(output.exists())
            self.assertEqual(subprocess.run(cmd+['--reference-2d', str(REFERENCE)], capture_output=True).returncode, 1)
            self.assertTrue(output.is_file())
            self.assertEqual(subprocess.run(cmd+['--reference-2d', str(REFERENCE)], capture_output=True).returncode, 2)

    def test_source_verification_rejects_unverified_archive_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)/'no.json'
            with self.assertRaisesRegex(ValueError, 'original fixture or archive bytes changed'):
                comparison.verify_source(REFERENCE, REFERENCE, REFERENCE, out)
            self.assertFalse(out.exists())


class ManufacturedFields(unittest.TestCase):
    def setUp(self):
        # Known affine field; these arrays are never written as reference data.
        x = np.linspace(0, 1, 661)
        self.target = np.column_stack((1+x, 2-x))

    def test_known_affine_pass_and_fail_on_fixed_five_percent_limit(self):
        for factor, expected in [(1., 'pass'), (1.049, 'pass'), (1.051, 'fail')]:
            result = comparison.displacement_metrics(self.target, {d: factor*self.target for d in comparison.DEPTHS})
            self.assertEqual(result['five_percent_verdict'], expected)
            for value in result['displacement_sample_relative_L2_per_depth'].values():
                self.assertAlmostEqual(value, abs(factor-1), places=14)
            self.assertEqual(result['displacement_sample_relative_L2_limit'], .05)

    def test_exact_five_percent_boundary_is_inclusive(self):
        target = np.zeros((661, 2)); target[0, 0] = 20.
        sample = target.copy(); sample[0, 0] = 21.
        result = comparison.displacement_metrics(target, {d: sample for d in comparison.DEPTHS})
        self.assertEqual(result['displacement_sample_relative_L2_per_depth']['0.5'], .05)
        self.assertEqual(result['five_percent_verdict'], 'pass')

    def test_front_or_back_failure_cannot_be_hidden_by_middle_plane(self):
        samples = {d: self.target.copy() for d in comparison.DEPTHS}
        samples['1.0'] *= 1.06
        result = comparison.displacement_metrics(self.target, samples)
        self.assertEqual(result['five_percent_verdict'], 'fail')
        self.assertEqual(result['displacement_sample_relative_L2_per_depth']['0.5'], 0.)
        self.assertFalse(result['depth_variation_passed'])

    def test_all_nodes_both_components_and_maximum_error_are_used(self):
        samples = {d: self.target.copy() for d in comparison.DEPTHS}
        for field in samples.values(): field[-1, 0] += 1.
        result = comparison.displacement_metrics(self.target, samples)
        self.assertEqual(result['five_percent_verdict'], 'pass')
        self.assertFalse(result['maximum_displacement_error_passed'])

    def test_missing_nonfinite_wrong_count_and_zero_norm_are_refused(self):
        valid = {d: self.target.copy() for d in comparison.DEPTHS}
        bad = []
        bad.append((np.zeros((661, 2)), valid))
        bad.append((self.target[:-1], valid))
        bad.append((self.target, {'0.5': self.target}))
        for value in (float('nan'), float('inf')):
            samples = copy.deepcopy(valid); samples['0.0'][0, 0] = value
            bad.append((self.target, samples))
        samples = copy.deepcopy(valid); samples['0.5'] = samples['0.5'][:-1]
        bad.append((self.target, samples))
        for target, samples in bad:
            with self.assertRaises(ValueError): comparison.displacement_metrics(target, samples)


class UserPcPreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.folder = Path(cls.temp.name)/'matched'
        cls.manifest = preparation.prepare(cls.folder, reference_path=REFERENCE)
        cls.sub = cls.folder/'xy1-ny1'
        cls.mesh = read(cls.sub/'mesh.json')
        cls.ground = read(cls.sub/'ground.json'); cls.case = read(cls.sub/'case.json')

    @classmethod
    def tearDownClass(cls): cls.temp.cleanup()

    def manufactured_run(self):
        state = dict(displacement_m=np.zeros((len(self.mesh['nodes']), 3)).tolist(),
                     active_elements=[True]*len(self.mesh['tetrahedra']))
        return dict(complete=True, ground=self.ground, case=self.case,
                    provenance=dict(solver_binary_sha256='synthetic-unit-test'),
                    geostatic=dict(complete=True, progress=1., last_committed=state, attempts=[],
                                   final_relative_physical_force_residual=0., maximum_absolute_engineering_strain_from_seed=0.),
                    frames=[dict(state=state, applied_force_kn=[0., 0., -4400.], support_reaction_kn=[0., 0., 4400.],
                                 relative_force_residual=0., maximum_normalized_yield=0.,
                                 maximum_absolute_engineering_strain_from_seed=0.), {}])

    def test_matching_material_exact_661_maps_and_unmodified_f1_defaults(self):
        self.assertFalse(self.manifest['solved'])
        self.assertFalse(self.manifest['griffiths_slope_five_percent_closed'])
        self.assertEqual([v['elements'] for v in self.manifest['levels']], [1200, 4800, 19200])
        ref = comparison.load_reference(REFERENCE)
        for level in self.manifest['levels']:
            sub = self.folder/level['folder']; mesh = read(sub/'mesh.json')
            xyz = np.asarray(mesh['nodes']); cells = np.asarray(mesh['tetrahedra'])
            self.assertAlmostEqual(float(np.linalg.det(np.transpose(xyz[cells[:, 1:]]-xyz[cells[:, :1]], (0, 2, 1))).sum()/6), 220.)
            mapping = read(sub/'reference-map.json')
            self.assertEqual(mapping, comparison.sample_nodes(mesh, ref))
            for depth, indices in mapping.items():
                self.assertEqual(len(indices), 661)
                self.assertTrue(np.allclose(xyz[indices][:, [0, 2]], np.asarray(ref['nodes'])[:, 1:3], rtol=0, atol=1e-12))
                self.assertTrue(np.all(xyz[indices][:, 1] == float(depth)))
            self.assertEqual(read(sub/'case.json'), read(sub/'stages.json'))
            for name, h in level['files_sha256'].items(): self.assertEqual(sha(sub/name), h)
        self.assertEqual(self.ground['materials'][0]['mohr_coulomb']['cohesion_kpa'], 7.407407407407406)
        self.assertEqual(self.ground['materials'][0]['mohr_coulomb']['friction_deg'], 15.088616807815535)
        baseline = ROOT/'interim/2026-11-07-plastic-spatial-convergence'
        self.assertTrue(baseline.is_dir())

    def test_missing_reference_fails_before_creating_output(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)/'never-created'
            with self.assertRaises(FileNotFoundError):
                preparation.prepare(out, reference_path=Path(temp)/'missing.json')
            self.assertFalse(out.exists())

    def test_no_interpolation_when_original_sample_node_missing(self):
        mesh = copy.deepcopy(self.mesh); mesh['nodes'][read(self.sub/'reference-map.json')['0.0'][0]][0] += .001
        with self.assertRaisesRegex(ValueError, 'missing original 2D node'):
            comparison.sample_nodes(mesh, read(REFERENCE))

    def test_incomplete_result_is_rejected_before_independent_audit(self):
        output = self.folder/'never-written.json'
        with patch.object(comparison, 'independent_audit') as audit:
            with self.assertRaises(FileNotFoundError):
                comparison.compare_series(self.folder, REFERENCE, output, REFERENCE)
            audit.assert_not_called()
        self.assertFalse(output.exists())

    def test_rehashed_manifest_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp); manifest = copy.deepcopy(self.manifest)
            manifest['levels'][1]['relative_h_xz'] = .3
            write(folder/'preparation.json', manifest)
            with self.assertRaisesRegex(ValueError, 'preparation manifest changed'):
                comparison.compare_series(folder, REFERENCE, folder/'no.json', REFERENCE)

    def test_f1_series_cannot_use_the_reduced_reference(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp); write(folder/'preparation.json', dict(material_reduction_factor=1.))
            with self.assertRaisesRegex(ValueError, 'F=1 cannot be compared'):
                comparison.compare_series(folder, REFERENCE, folder/'no.json', REFERENCE)

    def test_prepared_file_changes_cannot_be_hidden_by_original_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp); write(folder/'preparation.json', self.manifest)
            sub = folder/'xy1-ny1'; sub.mkdir(); write(sub/'mesh.json', {})
            with self.assertRaisesRegex(ValueError, 'prepared input changed'):
                comparison.compare_series(folder, REFERENCE, folder/'no.json', REFERENCE)

    def test_series_orchestration_requires_all_independent_audits(self):
        # Manufactured zero fields exercise orchestration only. No FEM solve,
        # and no synthetic result is retained outside this temporary directory.
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)/'series'; manifest = preparation.prepare(folder, reference_path=REFERENCE)
            for level in manifest['levels']:
                sub = folder/level['folder']; mesh = read(sub/'mesh.json')
                run = self.manufactured_run()
                state = dict(displacement_m=np.zeros((len(mesh['nodes']), 3)).tolist(),
                             active_elements=[True]*len(mesh['tetrahedra']))
                run['geostatic']['last_committed'] = state; run['frames'][0]['state'] = state
                run['provenance']['solver_binary_sha256'] = sha(REFERENCE)
                write(sub/'result.json', run)
            with patch.object(comparison, 'independent_audit', side_effect=AssertionError('audit failed')) as audit:
                with self.assertRaisesRegex(AssertionError, 'audit failed'):
                    comparison.compare_series(folder, REFERENCE, folder/'no.json', REFERENCE)
                self.assertEqual(audit.call_count, 1)
                self.assertFalse((folder/'no.json').exists())
            with patch.object(comparison, 'independent_audit') as audit:
                report = comparison.compare_series(folder, REFERENCE, folder/'synthetic-test-only.json', REFERENCE)
                self.assertEqual(audit.call_count, 3)
            self.assertEqual([v['five_percent_verdict'] for v in report['levels']], ['fail']*3)
            self.assertFalse(report['all_levels_five_percent_passed'])
            self.assertFalse(report['griffiths_slope_five_percent_closed'])
            self.assertEqual(report['continuum_verdict'], 'indeterminate')

    def test_completed_run_preflight_and_physical_rejections(self):
        run = self.manufactured_run()
        self.assertEqual(comparison.validate_run(run, self.mesh, self.ground, self.case, 'synthetic-unit-test').shape,
                         (len(self.mesh['nodes']), 3))
        for change in ('incomplete', 'force', 'global_force', 'strain', 'yield', 'transverse', 'nonfinite',
                       'wrong_binary', 'wrong_material', 'wrong_frame', 'inactive', 'rejected_history'):
            bad = copy.deepcopy(run)
            if change == 'incomplete': bad['geostatic']['complete'] = False
            if change == 'force': bad['frames'][0]['relative_force_residual'] = 1e-7
            if change == 'global_force': bad['frames'][0]['support_reaction_kn'][2] += 1.
            if change == 'strain': bad['frames'][0]['maximum_absolute_engineering_strain_from_seed'] = .051
            if change == 'yield': bad['frames'][0]['maximum_normalized_yield'] = 1e-6
            if change == 'transverse': bad['geostatic']['last_committed']['displacement_m'][0][1] = 1e-20
            if change == 'nonfinite': bad['geostatic']['last_committed']['displacement_m'][0][0] = float('nan')
            if change == 'wrong_binary': bad['provenance']['solver_binary_sha256'] = 'changed'
            if change == 'wrong_material': bad['ground']['materials'][0]['mohr_coulomb']['cohesion_kpa'] = 10.
            if change == 'wrong_frame': bad['frames'][0]['state'] = {}
            if change == 'inactive': bad['geostatic']['last_committed']['active_elements'][0] = False
            if change == 'rejected_history': bad['geostatic']['attempts'] = [dict(accepted=False, committed_state_unchanged_on_rejection=False)]
            with self.subTest(change=change), self.assertRaises(ValueError):
                comparison.validate_run(bad, self.mesh, self.ground, self.case, 'synthetic-unit-test')


if __name__ == '__main__': unittest.main()
