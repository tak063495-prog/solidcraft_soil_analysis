"""Light A-03 contracts only; full equilibrium evidence is sealed outside CI."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-12-22-a03-wet-soil-birth-cause-diagnosis'
sys.path.insert(0,str(HERE))
from a03cause_metrics import dot,norm,metric,compare,common
from a03cause_reference import volume_quadrature
from a03cause_driver import reverse_keys
from a03cause_binding import bind,study
import fem_reference as fem


class A03WetSoilBirthCauseDiagnosis(unittest.TestCase):
    def test_fixed_order_norm_keeps_small_denominator_observable(self):
        self.assertEqual(dot([1e16,1.,-1e16],[1.,1.,1.]),1.)
        self.assertEqual(norm([[3.,4.],[0.,12.]]),13.)
        value=metric([1e-18],[0.]);self.assertEqual(value['denominator'],1e-12)
        self.assertTrue(value['denominator_floor_active']);self.assertEqual(value['reference_l2'],0.)
        with self.assertRaises(ValueError):metric([1.],[1.,2.])
        with self.assertRaises(ValueError):metric([math.nan],[0.])

    def test_common_grid_does_not_interpolate(self):
        row=lambda p:(dict(stage='birth',progress=p),dict(u=np.array([p])))
        a=[row(.5),row(1.)];b=[row(.25),row(.5),row(.75),row(1.)]
        self.assertEqual(compare(a,common(b,2),set())['fields']['u']['absolute_l2'],0.)
        with self.assertRaises(ValueError):compare(a,b,set())
        with self.assertRaises(ValueError):compare(a,[row(.25),row(1.)],set())
        with self.assertRaises(ValueError):compare(a,[(a[0][0],{}),a[1]],set())

    def test_generic_import_names_are_restored(self):
        names=('fields','reference','study');before={n:sys.modules.get(n) for n in names}
        loaded=bind();self.assertEqual(set(loaded),set(names))
        for name in names:self.assertIs(sys.modules.get(name),before[name])

    def test_volume_control_restores_default_even_after_failure(self):
        original=fem.volume;mesh,_,_=study.fixture(2)
        with volume_quadrature(5):self.assertIs(fem.volume,original)
        with self.assertRaisesRegex(RuntimeError,'retained failure'):
            with volume_quadrature(7):
                self.assertEqual(len(fem.volume(mesh,0,independent=True)['weight']),343)
                self.assertEqual(len(fem.volume(mesh,0)['weight']),27)
                raise RuntimeError('retained failure')
        self.assertIs(fem.volume,original)
        self.assertEqual(len(fem.volume(mesh,0,independent=True)['weight']),125)
        for bad in (True,3,9):
            with self.assertRaises(ValueError):
                with volume_quadrature(bad):pass

    def test_object_key_order_is_not_a_physical_input_change(self):
        source={'level':1.5,'fresh':{'stress':[5.,5.,5.,0.,0.,0.],'history':0.}}
        changed=reverse_keys((source,))
        self.assertEqual(changed,[source]);self.assertEqual(list(changed[0]),['fresh','level'])
        self.assertEqual(list(changed[0]['fresh']),['history','stress'])
        self.assertEqual(list(source),['level','fresh'])

    def test_GP_only_water_test_cannot_authorize_partial_birth(self):
        mesh,_,_=study.fixture(2);zmax=max(p[2] for p in mesh['nodes']);level=math.nextafter(zmax,-math.inf)
        z=np.concatenate([fem.volume(mesh,e)['position'][:,2] for e in range(len(mesh['tetrahedra']))])
        self.assertLess(max(z),level);self.assertGreater(zmax,level)
        self.assertNotEqual(level,zmax)

    def test_document_terminal_status_with_real_crlf(self):
        text=(HERE/'README.ja.md').read_text(encoding='utf-8')
        self.assertTrue(text.endswith('A-03は未達、passed=false。\n'))
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'readme';p.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
            self.assertEqual(p.read_text(encoding='utf-8'),text)

    def test_sealed_parts_are_available_and_sha_fixed(self):
        a=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
        path=HERE/'archive-index.json';self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),a['archive_index_sha256'])
        data=json.loads(path.read_text(encoding='utf-8'))
        for part in data['parts']:
            p=HERE/part['name'];self.assertEqual(p.stat().st_size,part['size'])
            with p.open('rb') as f:self.assertEqual(hashlib.file_digest(f,'sha256').hexdigest(),part['sha256'])
        for required in ('formal/h1_5-n4/native-fields.jsonl','formal/h2_0-n4-gl7/independent-fields.jsonl','formal/refusals/darcy-before-birth-soil/flow-complete-checkpoint.json'):
            self.assertIn(required,data['entries'])

    def test_hydrostatic_accounting_stays_null_and_unachieved(self):
        a=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
        self.assertIs(a['passed'],False);self.assertIs(a['a03_complete'],False)
        data=json.loads((HERE/'measurements.json').read_text(encoding='utf-8'))
        self.assertIs(data['all_field_numerical_acceptance'],False)
        for key in ('flux','storage','birth_water_transfer','balance_acceptance'):self.assertIsNone(data['water_accounting'][key])
        self.assertEqual(len(data['refusals']),24);self.assertEqual(data['completed_flow_prefixes'],3)
        self.assertEqual(len(data['additional_refusals']),1)
        self.assertTrue(all(r['process']['exit_code']==1 and not r['result_file_exists'] and not r['process']['timed_out'] for r in data['refusals']))

    def test_saved_force_and_initialization_observations_remain_limited(self):
        data=json.loads((HERE/'measurements.json').read_text(encoding='utf-8'))
        for name,record in data['initialization'].items():
            self.assertTrue(record['fresh_state_equals_input_exact']);self.assertTrue(record['reordered_keys_preparation_and_final_state_exact'])
            self.assertTrue(record['repeated_native_run_exact_except_elapsed'])
            self.assertEqual(record['seeded_effective_stress'],[5.,5.,5.,0.,0.,0.])
        wet=data['force_paths']['h2_0-n4'];dry=data['force_paths']['h1_5-n4']
        self.assertLess(wet['maximum_buoyancy_identity_defect'],1e-10)
        self.assertGreater(dry['maximum_buoyancy_identity_defect'],1.)
        self.assertLess(dry['maximum_defect_plus_omitted_exterior'],1e-10)
        self.assertTrue(dry['omitted_exterior_is_force_accounting_not_a_new_solved_case'])
