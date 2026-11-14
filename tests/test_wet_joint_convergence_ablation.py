"""Actual saved histories must not be promoted by sampling, interpolation or software success."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-11-15-wet-joint-convergence'
sys.path.insert(0,str(HERE))
spec=importlib.util.spec_from_file_location('wet_convergence_evidence',HERE/'evidence.py')
evidence=importlib.util.module_from_spec(spec);spec.loader.exec_module(evidence)
import study as ablation
import measure as measurements


class WetJointConvergenceAblation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.addClassCleanup(cls.temp.cleanup)
        cls.folder=Path(cls.temp.name)/'archive';cls.manifest=evidence.extract(cls.folder)
        cls.report=evidence.validate_report(ablation.read(cls.folder/'measurements.json'))
        cls.paths=cls.folder/'paths'

    def test_global_parameter_and_nonmonotone_failure_are_preserved(self):
        self.assertIn('not mesh',self.report['original_16_128_parameter'])
        rows=self.report['global_only_cut']
        for pair,value in [('16__32',.6762689785372742),('32__64',.0366224927491334),('64__128',.07425902835843268)]:
            self.assertAlmostEqual(rows[pair]['fields']['plastic_slip']['pooled_relative_l2'],value,places=13)
            self.assertFalse(rows[pair]['full_history_acceptance'])
        self.assertLess(self.report['richardson']['global_cut']['fields']['plastic_slip']['observed_order'],0)

    def test_saved_native_path_and_input_only_reference_are_separate(self):
        for route in ablation.read(self.paths/'routes.json'):
            if route['native'] and route['id'] in self.report['matched_native_vs_independent']:
                folder=self.paths/route['id'];verification=ablation.read(folder/'native-raw/verification.json')
                self.assertEqual(verification['accepted_states'],7*route['g'])
                self.assertTrue(verification['uninterrupted_run_exact_except_elapsed'])
                self.assertFalse(ablation.read(folder/'independent-result.json')['reference_uses_native_results'])
                with (folder/'native-raw/accepted.jsonl').open() as f:
                    self.assertEqual(sum(1 for _ in f),7*route['g'])

    def test_contact_history_physics_are_checked_at_every_accepted_state(self):
        for name in ('h0-g64-l8-cut','h0-g128-l8-cut'):
            folder=self.paths/name
            result=ablation.audit(folder/'independent',ablation.read(self.paths/'mesh-h0.json'),ablation.read(self.paths/'ground.json'),ablation.read(self.paths/'stages.json'),'cut',ablation.read(folder/'configuration.json')['g'])
            expected=self.report['path_audits'][name+'/independent']
            self.assertEqual(result['accepted_states'],expected['accepted_states'])
            self.assertTrue(result['all_saved_states_audited'])
            self.assertFalse(result['numerical_acceptance'])
            self.assertEqual(result['existing_force_audit_tolerance'],expected['existing_force_audit_tolerance'])
            self.assertLess(result['maximum_reassembled_force_residual'],expected['existing_force_audit_tolerance'])
            self.assertEqual(result['maximum_pressure_change_kpa'],0)

    def test_missing_interior_chunk_cannot_pass_integrity_or_history(self):
        folder=self.paths/'h0-g128-l8-cut/independent';manifest=ablation.read(folder/'manifest.json');path=folder/manifest['chunks'][1]['file'];raw=path.read_bytes()
        try:
            path.unlink()
            with self.assertRaises(FileNotFoundError):list(ablation.states(folder))
        finally:path.write_bytes(raw)

    def test_changed_contact_field_is_rejected_even_if_chunk_is_resealed(self):
        folder=self.paths/'h0-g64-l8-cut/independent';mp=folder/'manifest.json';mr=mp.read_bytes();manifest=json.loads(mr);p=folder/manifest['chunks'][0]['file'];raw=p.read_bytes()
        try:
            with np.load(p,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
            if manifest.get('encoding'):
                decoded=ablation.codec.decode(arrays['traction'],manifest['encoding']);decoded[0,0,2]+=1
                arrays['traction']=ablation.codec.encode(decoded)
            else:arrays['traction'][0,0,2]+=1
            with p.open('wb') as f:np.savez_compressed(f,**arrays)
            manifest['chunks'][0]['sha256']=ablation.sha(p);mp.write_text(json.dumps(manifest))
            with self.assertRaises(AssertionError):
                ablation.audit(folder,ablation.read(self.paths/'mesh-h0.json'),ablation.read(self.paths/'ground.json'),ablation.read(self.paths/'stages.json'),'cut',64)
        finally:p.write_bytes(raw);mp.write_bytes(mr)

    def test_new_spatial_and_event_metrics_cannot_override_old_acceptance(self):
        for key in ('passed','all_history_accepted','p6_complete'):
            forged=copy.deepcopy(self.report);forged[key]=True
            with self.assertRaisesRegex(ValueError,'acceptance flag'):evidence.validate_report(forged)
        for row in self.report['space_only'].values():
            self.assertFalse(row['pointwise_contact_fields_compared']);self.assertFalse(row['spatial_convergence_accepted'])
        for row in self.report['event_aligned_diagnostic'].values():
            self.assertFalse(row['accepted']);self.assertFalse(row['existing_metrics_replaced'])


if __name__=='__main__':unittest.main()
