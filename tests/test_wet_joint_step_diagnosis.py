"""Protect complete C-01 evidence, finite-reference limits and unchanged gates."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-11-18-wet-joint-step-convergence'
spec=importlib.util.spec_from_file_location('c01_evidence_tests',HERE/'evidence_steps.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)

class WetJointStepDiagnosis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.addClassCleanup(cls.temp.cleanup)
        cls.folder=Path(cls.temp.name)/'evidence';e.extract(cls.folder)
        cls.report=e.validate_report(e.read(cls.folder/'measurements.json'))

    def test_recomputed_nonmonotone_history_is_not_accepted(self):
        a=self.folder/'paths/g64/fixed8';b=self.folder/'paths/g128/fixed8'
        actual=e.m.old.compare(a,b)
        self.assertEqual(actual,self.report['pairs']['fixed8/64__128'])
        self.assertLess(self.report['observed_orders']['fixed8/32__64__128']['fields']['plastic_slip']['observed_order'],0)
        for key in ('passed','c02_c03_c04_accepted','p6_complete'):
            changed=copy.deepcopy(self.report);changed[key]=True
            with self.assertRaisesRegex(ValueError,'unsupported'):e.validate_report(changed)

    def test_all_entities_and_event_acceptance_are_retained(self):
        fields=e.read(self.folder/'all-entity-differences.json')['native/128__256']['all_entities']
        for key,n in [('u',54),('stress',324),('plastic_strain',324),('plastic_slip',32),('accumulated',32)]:
            self.assertEqual([p['index'] for p in fields[key]],list(range(n)))
        events=e.read(self.folder/'switch-events.json')['native/256']
        opening=[v for v in events['events'] if v['kind']=='opening' and v['stage']=='open']
        self.assertEqual({v['point'] for v in opening},set(range(32)))
        for v in events['events']:
            a=v['attempt'];self.assertTrue(a['accepted'])
            self.assertLessEqual(a['relative_force_residual'],a['force_limit'])
            self.assertLessEqual(a['maximum_true_linear_residual'],a['linear_limit'])

    def test_real_file_checkpoints_match_all_suffix_states(self):
        restart=self.report['separate_process_restart'];self.assertTrue(restart['all_seven_phases'])
        self.assertTrue(restart['uninterrupted_exact_except_elapsed'])
        direct=e.m.st.without_clock(e.read(self.folder/'paths/g256/native-raw/run.json'))
        for row in restart['checks']:
            target=self.folder/'paths/g256'/('restart-'+Path(row['checkpoint']).stem)
            self.assertEqual(e.m.st.without_clock(e.read(target/'run.json')),direct)
            self.assertGreater(row['suffix_states_exact'],0)
        self.assertEqual(len({row['process_id'] for row in restart['checks']}),7)

    def test_archived_sources_import_in_isolated_python(self):
        script=self.folder/'sources/interim/2026-11-18-wet-joint-step-convergence/measure_steps.py'
        result=subprocess.run([sys.executable,'-I','-X','utf8',str(script),'--help'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_missing_interior_field_cannot_pass_as_complete(self):
        folder=self.folder/'paths/g256/native';manifest=e.read(folder/'manifest.json');p=folder/manifest['chunks'][1]['file']
        raw=p.read_bytes()
        try:
            p.unlink()
            with self.assertRaises(FileNotFoundError):list(e.m.states(folder))
        finally:p.write_bytes(raw)

    def test_contact_replay_and_chord_defect_are_distinct_from_solution(self):
        for row in self.report['frozen_history_replay'].values():
            self.assertTrue(row['all_newton_iterates_share_immutable_old_history'])
            self.assertTrue(row['material_equivalent_plastic_strain_and_work_zero'])
        chord=e.read(self.folder/'chord-defect.json')
        self.assertFalse(chord['numerical_solution'] or chord['physical_coordinates_changed'] or chord['continuous_reference'])
        self.assertGreater(chord['stage_maxima']['open']['maximum_interior_force_residual'],1e-8)
        self.assertLess(chord['stage_maxima']['open']['maximum_endpoint_force_residual'],1e-10)
        self.assertEqual(chord['local_partitions_chord'],chord['local_partitions_curved'])

if __name__=='__main__':unittest.main()
