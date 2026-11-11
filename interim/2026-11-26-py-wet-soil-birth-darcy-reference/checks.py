"""Explicit small suite. Filename deliberately outside existing test_*.py discovery."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from reference import (solve, export, exact_surface, exact_head, conductivity, birth,
                       policies, GW, RECHARGE)
from compare import diagnose, read, read_json, write_csv


class FreeBoundary(unittest.TestCase):
    def test_unknown_surface_inside_new_body(self):
        for a in (.5, -.25):
            self.assertTrue(.7 < solve(16, a)['surface'] < 2.)

    def test_pde_closed_form_independent_derivative(self):
        for a in (.5, -.25):
            z = np.linspace(.1, .9, 9)
            derivative = (exact_head(z+1e-5, a)-exact_head(z-1e-5, a))/2e-5
            np.testing.assert_allclose(-conductivity(z, a)*derivative, -RECHARGE, rtol=1e-8)

    def test_surface_boundary(self):
        for a in (.5, -.25):
            s = exact_surface(a)
            self.assertAlmostEqual(float(exact_head(s, a)), s, places=12)
            r = solve(32, a)
            self.assertLess(abs(r['head'][-1]-r['surface']), 1e-11)
            self.assertTrue(np.all(r['head']-r['z'] >= -1e-11))

    def test_second_order_surface_not_roundoff(self):
        for a in (.5, -.25):
            errors = [abs(solve(n, a)['surface']-exact_surface(a)) for n in (8,16,32,64)]
            self.assertGreater(errors[-1], 1e-8)
            self.assertTrue(all(1.99 < np.log2(x/y) < 2.01 for x,y in zip(errors, errors[1:])))

    def test_head_convergence_separate_metric(self):
        for a in (.5, -.25):
            values = [export(n, a)['metrics']['head_rms_m'] for n in (8,16,32,64)]
            self.assertTrue(all(1.99 < np.log2(x/y) < 2.04 for x,y in zip(values,values[1:])))

    def test_local_and_global_conservation(self):
        for a in (.5, -.25):
            r = export(64,a)
            self.assertLess(abs(r['metrics']['water_defect_m3_s']), 1e-16)
            self.assertLess(r['metrics']['local_flux_spread_m3_s'], 1e-16)
            self.assertGreater(r['budget'][0]['outward_m3_s'], 0.)
            self.assertLess(r['budget'][1]['outward_m3_s'], 0.)

    def test_no_negative_pressure_clipping_in_dry_region(self):
        r = export(16)
        dry = [s for s in r['samples'] if s['z_m'] > r['surface_m']]
        self.assertTrue(dry)
        self.assertTrue(all(s['pressure_kpa'] == 0 for s in dry))
        self.assertFalse(r['metadata']['pressure_clipped'])

    def test_bad_grid(self):
        for n in (0, 3, True, 3.5, 20000):
            with self.assertRaises(ValueError): solve(n)

    def test_bad_material(self):
        with self.assertRaises(ValueError): solve(8, 0.)


class BirthBalance(unittest.TestCase):
    def test_all_node_equilibrium(self):
        for a in (.5, -.25):
            for n in (8, 16, 32, 64):
                self.assertLess(export(n,a)['metrics']['nodal_residual_max_kn'], 1e-11)

    def test_pore_divergence_identity(self):
        self.assertLess(export()['metrics']['pore_identity_max_kn'], 1e-11)

    def test_buoyancy_seepage_separate_not_double_counted(self):
        r=export()
        self.assertLess(r['metrics']['buoyancy_seepage_identity_max_kn'], 1e-11)
        self.assertGreater(sum(f['buoyancy_kn'] for f in r['forces']), 0)
        self.assertLess(sum(f['seepage_kn'] for f in r['forces']), 0)

    def test_independent_global_total_weight(self):
        r=export()
        self.assertAlmostEqual(sum(f['weight_kn'] for f in r['forces']),
                               -20*r['surface_m']-17*(2-r['surface_m']), places=11)

    def test_total_effective_pressure_identity(self):
        for s in export()['samples']:
            self.assertAlmostEqual(s['total_kpa']-s['effective_kpa'], s['pressure_kpa'], places=12)

    def test_regenerated_faces_no_stale_exterior(self):
        t=export()['metrics']['topology']
        self.assertEqual(t['stale_old_exterior_after_kn'],0)
        self.assertEqual(t['paired_interface_sum_kn'],0)
        self.assertEqual(t['before_old_exposed_kn'],t['after_old_interface_kn'])
        self.assertLess(t['after_old_interface_kn'],0)

    def test_policies_distinct_roles(self):
        p=policies(5.,1.2,1.2,1.2)
        self.assertEqual(p['soil']['normal_axis'],'x')
        self.assertEqual(p['soil']['structure_pore_pressure_kpa'],0)
        self.assertEqual(p['continuous']['normal_traction_kn'],[-5.,5.])
        self.assertEqual(p['prescribed']['fixed_head'],1.2)

    def test_policy_conflicts_not_averaged(self):
        with self.assertRaises(ValueError): policies(5.,1.2,1.3,1.2)


class Diagnostics(unittest.TestCase):
    def setUp(self):
        self.reference=export(8)
        self.data=copy.deepcopy(self.reference)

    def check(self): return diagnose(self.data,self.reference)

    def test_identity_not_acceptance(self):
        r=self.check()
        self.assertIsNone(r['accepted'])
        self.assertFalse(r['a03_complete'])

    def test_csv_roundtrip(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'csv'
            write_csv(self.data,path)
            self.assertEqual(diagnose(read(path),self.reference),self.check())

    def test_wrong_units(self):
        self.data['metadata']['units']['pressure']='Pa'
        with self.assertRaises(ValueError): self.check()

    def test_wrong_sign(self):
        self.data['metadata']['stress_sign']='tension_positive'
        with self.assertRaises(ValueError): self.check()

    def test_wrong_flow_sign(self):
        self.data['metadata']['flow_sign']='inward_positive'
        with self.assertRaises(ValueError): self.check()

    def test_wrong_force_sign(self):
        self.data['metadata']['force_sign']='compression_positive'
        with self.assertRaises(ValueError): self.check()

    def test_wrong_stress_type(self):
        self.data['metadata']['stress_kind']='total'
        with self.assertRaises(ValueError): self.check()

    def test_pressure_head_not_total_head(self):
        self.data['metadata']['head_kind']='pressure'
        with self.assertRaises(ValueError): self.check()

    def test_missing_id(self):
        self.data['samples'].pop()
        with self.assertRaises(ValueError): self.check()

    def test_duplicate_id(self):
        self.data['samples'].append(self.data['samples'][0])
        with self.assertRaises(ValueError): self.check()

    def test_extra_id(self):
        r=copy.deepcopy(self.data['samples'][0]); r['id']='extra'; self.data['samples'].append(r)
        with self.assertRaises(ValueError): self.check()

    def test_different_mesh_refused(self):
        self.data['samples'][0]['z_m']+=.001
        with self.assertRaises(ValueError): self.check()

    def test_interpolation_refused(self):
        self.data['metadata']['interpolation']=True
        with self.assertRaises(ValueError): self.check()

    def test_filled_pressure_refused(self):
        self.data['metadata']['pressure_filled']=True
        with self.assertRaises(ValueError): self.check()

    def test_clipped_pressure_refused(self):
        self.data['metadata']['pressure_clipped']=True
        with self.assertRaises(ValueError): self.check()

    def test_missing_pressure(self):
        del self.data['samples'][0]['pressure_kpa']
        with self.assertRaises(ValueError): self.check()

    def test_nan(self):
        self.data['samples'][0]['pressure_kpa']=float('nan')
        with self.assertRaises(ValueError): self.check()

    def test_bool_not_force(self):
        self.data['forces'][0]['weight_kn']=False
        with self.assertRaises(ValueError): self.check()

    def test_incomplete_not_success(self):
        for status in ('timed_out','cancelled','running'):
            self.data['metadata']['status']=status
            with self.assertRaises(ValueError): self.check()

    def test_balanced_force_corruption_visible_per_node(self):
        self.data['forces'][0]['weight_kn']+=1
        self.data['forces'][1]['weight_kn']-=1
        r=self.check()
        self.assertEqual(r['differences']['forces']['weight_kn']['max_abs'],1)
        self.assertGreater(r['nodal_equilibrium_max_kn'],.99)
        self.assertIsNone(r['accepted'])

    def test_false_numeric_label_identity_refused(self):
        self.data['samples'][0]['effective_kpa']=self.data['samples'][0]['total_kpa']
        with self.assertRaises(ValueError): self.check()

    def test_pressure_head_mislabeled_numerically_refused(self):
        self.data['hydraulic'][0]['head_m']-=.5
        with self.assertRaises(ValueError): self.check()

    def test_water_defect_not_hidden(self):
        self.data['budget'][0]['outward_m3_s']+=1e-5
        self.assertGreater(self.check()['outward_water_defect_m3_s'],9e-6)

    def test_duplicate_json_key(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'a.json'; p.write_text('{"x":0,"x":1}')
            with self.assertRaises(ValueError): read_json(p)

    def test_duplicate_csv_column(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'csv';write_csv(self.data,p)
            (p/'hydraulic.csv').write_text('id,id\na,b\n')
            with self.assertRaises(ValueError): read(p)


if __name__ == '__main__':
    unittest.main(verbosity=2)
