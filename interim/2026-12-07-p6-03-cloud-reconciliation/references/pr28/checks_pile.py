"""Explicit opt-in tests, outside existing tests/ and default test_*.py patterns."""
import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss

import pile_reference as ref
from compare_pile import ContractError, compare, index_rows, load, read_json, validate_metadata
from export_reference import export, verification_metadata


class ReferenceChecks(unittest.TestCase):
    def test_axial_end_conditions(self):
        for tip in [0., 0.7, 10.]:
            fields = ref.axial_exact([0., 2.], 2., 3., 4., 5., tip)
            np.testing.assert_allclose(fields[0, 1], 5., rtol=1e-13)
            np.testing.assert_allclose(fields[1, 1], tip*fields[1, 0], atol=1e-13)

    def test_axial_load_balance(self):
        g,w = leggauss(32)
        values = ref.axial_exact(g+1, 2., 3., 4., 5., 0.7)
        tip = ref.axial_exact([2.], 2., 3., 4., 5., 0.7)[0,1]
        self.assertAlmostEqual(float(-w@values[:,2]+tip), 5., places=12)

    def test_axial_known_exponential(self):
        s=np.linspace(0,2,9)
        f=ref.axial_exact(s,2.,3.,3.,5.,3.)
        np.testing.assert_allclose(f[:,0], 5/3*np.exp(-s), rtol=1e-13)
        np.testing.assert_allclose(f[:,1], 5*np.exp(-s), rtol=1e-13)

    def test_axial_orders(self):
        series=ref.convergence("axial",2.,3.,4.,5.,0.7)
        p=series[-1]["observed_order"]
        self.assertTrue(1.9 < p["u"] < 2.1)
        self.assertTrue(0.9 < p["N"] < 1.1)

    def test_beam_end_conditions(self):
        f=ref.beam_exact([0.,2.],2.,3.,4.,5.)
        np.testing.assert_allclose(f[:,2],[0.,0.],atol=1e-12)
        np.testing.assert_allclose(f[:,3],[-5.,0.],atol=1e-12)

    def test_beam_load_and_moment_balance(self):
        g,w=leggauss(32)
        s=g+1
        f=ref.beam_exact(s,2.,3.,4.,5.)
        self.assertAlmostEqual(float(-w@f[:,4]),5.,places=11)
        self.assertAlmostEqual(float(w@(s*f[:,4])),0.,places=11)

    def test_beam_orders(self):
        p=ref.convergence("beam",2.,3.,4.,5.)[-1]["observed_order"]
        for quantity, expected in [("w",4),("slope",3),("M",2),("V",1)]:
            self.assertLess(abs(p[quantity]-expected),0.15)

    def test_semi_infinite_end_and_decay(self):
        f=ref.beam_semi_infinite([0.,100.],3.,4.,5.)
        beta=(4/(4*3))**.25
        self.assertAlmostEqual(f[0,0],5/(2*3*beta**3))
        np.testing.assert_allclose(f[0,2:4],[0.,-5.],atol=1e-13)
        self.assertLess(np.max(np.abs(f[1])),1e-25)

    def test_finite_is_not_semi_infinite(self):
        finite=ref.beam_exact([0.],1.,1.,4.,1.)[0,0]
        infinite=ref.beam_semi_infinite([0.],1.,4.,1.)[0,0]
        self.assertGreater(abs(finite-infinite),.1)

    def test_fe_equilibrium(self):
        for kind in ["axial","beam"]:
            sol,residual=ref.finite_elements(kind,8,2.,3.,4.,5.)
            self.assertLess(residual,1e-10)
            g,w=leggauss(4)
            s=np.concatenate([(e+(g+1)/2)/4 for e in range(8)])
            weights=np.tile(w/8,8)
            reaction=ref.evaluate(kind,sol,s,2.,3.,4.)[:,-1]
            self.assertAlmostEqual(float(-weights@reaction),5.,places=9)

    def test_input_guards(self):
        for v in [0.,-1.,float('nan'),float('inf')]:
            with self.subTest(v=v), self.assertRaises(ValueError):
                ref.axial_exact([0.],1.,1.,v,1.,1.)
        with self.assertRaises(ValueError):
            ref.beam_exact([2.],1.,1.,4.,1.)
        with self.assertRaises(ValueError):
            ref.finite_elements("beam",8,1.,1.,4.,1.,1.)
        with self.assertRaises(ValueError):
            ref.finite_elements("axial",True,1.,1.,4.,1.)

    def test_undefined_acceptance(self):
        report=ref.report()
        self.assertIsNone(report['accepted'])
        self.assertFalse(report['p6_complete'])
        self.assertEqual(report['original_3d_comparison'],'判定不能')


class ComparisonChecks(unittest.TestCase):
    def setUp(self):
        self.meta=verification_metadata("axial")
        self.reference=export(self.meta,"axial","exact")
        self.actual=copy.deepcopy(self.reference)

    def result(self):
        return compare((self.reference['metadata'],self.reference['rows']),
                       (self.actual['metadata'],self.actual['rows']))

    def test_identical_has_no_acceptance(self):
        r=self.result()
        self.assertIsNone(r['accepted'])
        self.assertEqual(r['quantities']['N']['max_absolute_error'],0.)

    def test_row_order_irrelevant(self):
        self.actual['rows'].reverse()
        self.assertEqual(self.result()['quantities']['u']['relative_l2'],0.)

    def test_known_difference(self):
        for row in self.actual['rows']:
            row['value'] *= 1.1
        self.assertAlmostEqual(self.result()['quantities']['u']['relative_l2'],.1)

    def test_missing_id_rejected(self):
        self.actual['rows'].pop()
        with self.assertRaises(ContractError): self.result()

    def test_duplicate_id_rejected(self):
        self.actual['rows'].append(self.actual['rows'][0])
        with self.assertRaises(ContractError): self.result()

    def test_unknown_id_rejected(self):
        self.actual['rows'][0]['id']='other'
        with self.assertRaises(ContractError): self.result()

    def test_units_rejected(self):
        self.actual['metadata']['units']='mm-N-MPa'
        with self.assertRaises(ContractError): self.result()

    def test_quantity_units_rejected(self):
        self.actual['metadata']['quantities']['N']['unit']='kPa'
        with self.assertRaises(ContractError): self.result()

    def test_sign_rejected(self):
        self.actual['metadata']['signs']['N']='tension_positive'
        with self.assertRaises(ContractError): self.result()

    def test_basis_rejected(self):
        self.actual['metadata']['basis']['axes']=[[0,1,0],[-1,0,0],[0,0,1]]
        with self.assertRaises(ContractError): self.result()

    def test_left_handed_basis_rejected(self):
        self.actual['metadata']['basis']['axes'][2][2]=-1
        with self.assertRaises(ContractError): self.result()

    def test_coordinate_mismatch_rejected(self):
        self.actual['rows'][0]['z_m']+=.01
        with self.assertRaises(ContractError): self.result()

    def test_changed_physics_rejected(self):
        for key in ['geometry','materials','loads','boundary','interface','water','reference_state','model']:
            with self.subTest(key=key):
                self.actual=copy.deepcopy(self.reference)
                self.actual['metadata']['case'][key]='changed'
                with self.assertRaises(ContractError): self.result()

    def test_unknown_physics_rejected_even_on_both_sides(self):
        for value in [None,'未確認',{},'unknown']:
            with self.subTest(value=value):
                for envelope in [self.reference,self.actual]:
                    envelope['metadata']['case']['boundary']=value
                with self.assertRaises(ContractError): self.result()

    def test_changed_frame_rejected(self):
        self.actual['metadata']['frame']='total_gravity_plus_load'
        with self.assertRaises(ContractError): self.result()

    def test_nonfinite_rejected(self):
        for value in [float('nan'),float('inf'),'NaN','']:
            self.actual['rows'][0]['value']=value
            with self.subTest(value=value),self.assertRaises(ContractError): self.result()

    def test_empty_rows_rejected(self):
        self.actual['rows']=[]
        with self.assertRaises(ContractError): self.result()

    def test_csv_and_json_roundtrip(self):
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory)
            (d/'ref.json').write_text(json.dumps(self.reference),encoding='utf-8')
            (d/'meta.json').write_text(json.dumps(self.meta),encoding='utf-8')
            with (d/'out.csv').open('w',newline='',encoding='utf-8') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(self.reference['rows'][0]))
                writer.writeheader();writer.writerows(self.reference['rows'])
            r=compare(load(d/'ref.json'),load(d/'out.csv',d/'meta.json'))
            self.assertEqual(r['quantities']['N']['count'],17)
            with self.assertRaises(ContractError): load(d/'out.csv')

    def test_duplicate_json_keys_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'bad.json'
            path.write_text('{"case":1,"case":2}',encoding='utf-8')
            with self.assertRaises(ContractError):read_json(path)

    def test_reference_export_refuses_3d(self):
        self.meta['case']['model']='3d-solid-elastic-bonded'
        with self.assertRaises(ContractError):export(self.meta,'axial','exact')

    def test_beam_export_and_fe_comparison(self):
        meta=verification_metadata('beam')
        exact=export(meta,'beam','exact')
        fe=export(meta,'beam','fe',16)
        r=compare((exact['metadata'],exact['rows']),(fe['metadata'],fe['rows']))
        self.assertLess(r['quantities']['w']['relative_l2'],1e-6)
        self.assertIsNone(r['accepted'])

    def test_quantity_key_order_does_not_change_export(self):
        expected=export(self.meta,'axial','exact')
        self.meta['quantities']=dict(reversed(list(self.meta['quantities'].items())))
        actual=export(self.meta,'axial','exact')
        self.assertEqual(expected['rows'],actual['rows'])

    def test_missing_quantity_inventory_rejected(self):
        del self.actual['metadata']['expected_ids']['N']
        with self.assertRaises(ContractError):self.result()

    def test_missing_physical_condition_rejected(self):
        del self.actual['metadata']['case']['water']
        with self.assertRaises(ContractError):self.result()


class SolidChecks(unittest.TestCase):
    def test_hex8_linear_patch(self):
        import solid_pile_reference as solid
        widths=np.array([2.,3.,4.])
        xyz=(solid.SIGNS+1)*widths/2
        gradient=np.array([[.1,.2,.3],[.4,.5,.6],[.7,.8,.9]])
        u=(xyz@gradient.T+np.array([2.,3.,4.])).ravel()
        expected=np.array([.1,.5,.9,.6,1.4,1.])
        for point in [[0.,0.,0.],[.3,-.6,.7]]:
            n,b=solid.shape(point,widths)
            self.assertAlmostEqual(float(n.sum()),1.)
            np.testing.assert_allclose(b@u,expected,atol=1e-14)

    def test_hex8_rigid_motion_and_energy(self):
        import solid_pile_reference as solid
        widths=np.array([2.,3.,4.]);d=solid.tensor(100.,.25)
        xyz=(solid.SIGNS+1)*widths/2
        k=solid.element(widths,d)
        rotation=np.cross(np.array([.2,.3,.4]),xyz)+1
        np.testing.assert_allclose(k@rotation.ravel(),0.,atol=1e-11)
        u=xyz.ravel()*.1
        strain=np.array([.1,.1,.1,0.,0.,0.])
        self.assertAlmostEqual(float(u@k@u),float(strain@d@strain*np.prod(widths)),places=10)
        np.testing.assert_allclose(k,k.T,atol=1e-12)

    def test_same_3d_physics_and_balance(self):
        import solid_pile_reference as solid
        r=solid.solve(1)
        np.testing.assert_allclose(r['volumes_m3'],[794.,6.],atol=1e-12)
        self.assertEqual(r['head_area_m2'],1.)
        np.testing.assert_allclose(r['applied_gravity_kn'],[0.,0.,-16024.],atol=1e-9)
        np.testing.assert_allclose(r['applied_increment_kn'],[10.,0.,-100.],atol=1e-12)
        applied=np.column_stack((r['applied_gravity_kn'],r['applied_increment_kn']))
        np.testing.assert_allclose(np.array(r['support_reaction_kn'])+applied,0.,atol=1e-8)
        self.assertLess(max(r['free_relative_residual']),1e-10)
        self.assertIsNone(r['accepted'])

    def test_mesh_nested_and_pile_geometry(self):
        import solid_pile_reference as solid
        coarse=solid.mesh(1);fine=solid.mesh(2)
        self.assertEqual(len(fine[1]),8*len(coarse[1]))
        self.assertTrue(set(map(tuple,coarse[0])).issubset(set(map(tuple,fine[0]))))
        for xyz,cells,pile,widths in [coarse,fine]:
            centroids=xyz[cells].mean(axis=1)
            expected=((centroids[:,0]>4.5)&(centroids[:,0]<5.5)&
                      (centroids[:,1]>4.5)&(centroids[:,1]<5.5)&(centroids[:,2]>2.))
            np.testing.assert_array_equal(pile,expected)

    def test_3d_refuses_nonuniform_refinement_sequence(self):
        import solid_pile_reference as solid
        with self.assertRaises(ValueError):solid.report((1,2,3))


class EvidenceChecks(unittest.TestCase):
    def test_historical_values_and_judgement_preserved(self):
        import diagnose_records
        root=Path(__file__).resolve().parents[2]
        r=diagnose_records.diagnose(root)
        self.assertEqual(r['historical_tet10_counts'],[436,819,2279])
        self.assertEqual(r['historical_settlement_mm'],
                         [0.32521742815831706,0.3485245109509388,0.36577978075570766])
        np.testing.assert_allclose(r['recomputed_relative_difference'],
                                   [.06687358294838157,.04717393008743999],rtol=0.,atol=1e-14)
        self.assertFalse(r['historical_spatial_convergence_accepted'])
        self.assertIsNone(r['accepted'])

    def test_changed_historical_source_is_refused(self):
        import diagnose_records
        root=Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as directory:
            temporary=Path(directory)
            for name in diagnose_records.EVIDENCE:
                target=temporary/name;target.parent.mkdir(parents=True,exist_ok=True)
                target.write_bytes((root/name).read_bytes())
            target=temporary/diagnose_records.EVIDENCE[0]
            target.write_bytes(target.read_bytes()+b'\nchanged')
            with self.assertRaises(ValueError):diagnose_records.diagnose(temporary)

    def test_no_project_solver_imports_or_subprocesses(self):
        import ast
        folder=Path(__file__).resolve().parent
        allowed={'argparse','json','pathlib','time','hashlib','csv','math','re','numpy','scipy',
                 'pile_reference','compare_pile','export_reference'}
        for name in ['pile_reference.py','solid_pile_reference.py','compare_pile.py',
                     'export_reference.py','diagnose_records.py']:
            tree=ast.parse((folder/name).read_text(encoding='utf-8'))
            for node in ast.walk(tree):
                if isinstance(node,ast.Import):
                    modules=[alias.name.split('.')[0] for alias in node.names]
                elif isinstance(node,ast.ImportFrom):
                    modules=[node.module.split('.')[0]]
                else:continue
                self.assertTrue(set(modules)<=allowed,(name,modules))


if __name__=='__main__':
    unittest.main()
