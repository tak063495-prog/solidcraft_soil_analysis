"""Native schema roundtrip, strict FreshNative guards, and linear diagnostics."""
from copy import deepcopy
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))


class PcRerunFix(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(['cargo', 'build', '--locked', '--release', '--example', 'nc_input_roundtrip',
                        '--example', 'nc_dense_lu_diagnostic'], cwd=ROOT, check=True)
        suffix = '.exe' if sys.platform == 'win32' else ''
        cls.roundtrip = ROOT/f'target/release/examples/nc_input_roundtrip{suffix}'
        cls.dense = ROOT/f'target/release/examples/nc_dense_lu_diagnostic{suffix}'

    def inputs(self, denom=128):
        from probe_nc_midpoint import inputs
        values = list(inputs('midpoint'))
        values[2]['stages'][2]['initial_increment'] = 1/denom
        return values

    def normalized(self, values):
        child = subprocess.run([str(self.roundtrip)], input=json.dumps(values),
                               text=True, capture_output=True, check=True)
        return json.loads(child.stdout)

    def test_native_schema_roundtrip_and_both_original_increments(self):
        from nc_fixed_inputs import canonical_inputs, fixed_inputs_match
        for denom in (128, 256):
            with self.subTest(denom=denom):
                raw = self.inputs(denom); before = deepcopy(raw)
                normalized = self.normalized(raw)
                self.assertNotEqual(raw, normalized)
                self.assertTrue(fixed_inputs_match(normalized, raw))
                self.assertEqual(raw, before)
                canonical = canonical_inputs(*raw)
                self.assertEqual(canonical, canonical_inputs(*canonical))

    def test_history_input_guard_accepts_native_roundtrip_before_state_guard(self):
        from nc_fold_fresh import FreshNative
        with tempfile.TemporaryDirectory() as directory:
            for denom in (128, 256):
                folder = Path(directory)/str(denom); folder.mkdir()
                mesh, ground, case = self.normalized(self.inputs(denom))
                state = dict(stage_index=2, stage_progress=.95, active_elements=[])
                (folder/'mesh.json').write_text(json.dumps(mesh), encoding='utf-8')
                (folder/'run.json').write_text(json.dumps(dict(ground=ground, case=case, last_committed=state)), encoding='utf-8')
                (folder/'accepted.jsonl').write_text('\n'.join(json.dumps(dict(state=state)) for _ in range(2)), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'no accepted point'):
                    FreshNative(folder, Path('must-not-launch'), folder/'material')

    def test_all_physical_and_control_edits_remain_rejected(self):
        from nc_fixed_inputs import fixed_inputs_match
        original = self.inputs()
        edits = [
            ('geometry', lambda m,g,c: m['nodes'][0].__setitem__(0, .01)),
            ('material', lambda m,g,c: g['materials'][0]['mcc'].__setitem__('kappa', .021)),
            ('constraint gradient', lambda m,g,c: c['boundary']['constraints'][0].__setitem__('gradient', [1e-30,0.,0.])),
            ('traction gradient', lambda m,g,c: c['boundary']['surface_loads'][0].__setitem__('gradient_kpa_per_m', [[1e-30,0.,0.],[0.,0.,0.],[0.,0.,0.]])),
            ('pressure gradient', lambda m,g,c: g['initialization']['geostatic']['base_loads']['surface_loads'][0].__setitem__('gradient_kpa_per_m', [0.,0.,1e-30])),
            ('linear tolerance', lambda m,g,c: c['settings'].__setitem__('relative_linear_tolerance', math.nextafter(1e-10, math.inf))),
            ('force tolerance', lambda m,g,c: c['settings'].__setitem__('relative_force_tolerance', 1e-7)),
            ('solver', lambda m,g,c: c['settings'].__setitem__('linear_solver', {'method':'sparse_gmres_ilu0'})),
            ('minimum increment', lambda m,g,c: c['settings'].__setitem__('minimum_increment', 1e-7)),
            ('Newton limit', lambda m,g,c: c['settings'].__setitem__('maximum_newton_iterations', 61)),
            ('integration', lambda m,g,c: c['settings'].__setitem__('mcc_integration', 'backward_euler')),
        ]
        for component in range(3):
            edits.append((f'holding load component {component}', lambda m,g,c,j=component: g['initialization']['geostatic']['base_loads']['surface_loads'][1]['traction_kpa'].__setitem__(j, .125)))
        from nc_fold_fresh import FreshNative
        with tempfile.TemporaryDirectory() as directory:
            for index, (name, edit) in enumerate(edits):
                with self.subTest(edit=name):
                    changed = self.normalized(original); edit(*changed)
                    self.assertFalse(fixed_inputs_match(changed, original))
                    folder = Path(directory)/str(index); folder.mkdir()
                    mesh, ground, case = changed
                    (folder/'mesh.json').write_text(json.dumps(mesh), encoding='utf-8')
                    (folder/'run.json').write_text(json.dumps(dict(ground=ground, case=case)), encoding='utf-8')
                    with self.assertRaisesRegex(ValueError, 'fixed inputs'):
                        FreshNative(folder, Path('must-not-launch'), folder/'material')

    def test_unknown_fields_and_false_numeric_defaults_are_not_ignored(self):
        from nc_fixed_inputs import fixed_inputs_match
        original = self.inputs()
        for path, value in [('unknown', 0.), ('gradient', [False, 0., 0.])]:
            changed = self.normalized(original)
            changed[2]['boundary']['constraints'][0][path] = value
            self.assertFalse(fixed_inputs_match(changed, original))
        for value in (float('nan'), float('inf')):
            changed = self.normalized(original)
            changed[2]['settings']['relative_linear_tolerance'] = value
            self.assertFalse(fixed_inputs_match(changed, original))

    def test_decimal_residual_detects_cancellation(self):
        from diagnose_nc_dense_lu import decimal_residual
        residual = decimal_residual([[1e16,1.,-1e16],[0.,1.,0.],[0.,0.,1.]], [0.,1.,1.], [1.,1.,1.])
        self.assertAlmostEqual(residual, 1/math.sqrt(2), places=15)

    def test_offline_refinement_records_every_step_and_preserves_source(self):
        from diagnose_nc_dense_lu import run
        system = dict(matrix=[[2.,3.],[0.,4.]], rhs=[8.,8.], solution=[1.,2.], native_true_relative_residual=0.)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); path = root/'systems.jsonl'
            path.write_text(json.dumps(system)+'\n', encoding='utf-8'); before = path.read_bytes()
            report = run(path, self.dense, root/'audit', refinements=3)
            self.assertEqual(path.read_bytes(), before)
            self.assertEqual(report['tolerance'], 1e-10)
            self.assertFalse(report['numerical_acceptance'])
            self.assertFalse(report['refinement_applied_to_production'])
            rows = report['systems'][0]['observations']
            self.assertEqual([r['refinement'] for r in rows], [0,1,2,3])
            self.assertTrue(all(r['solution']==[1.,2.] and r['decimal80_relative_residual']==0. for r in rows))



if __name__ == '__main__':
    unittest.main()
