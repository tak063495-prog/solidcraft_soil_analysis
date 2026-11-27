"""Explicitly invoked unittest suite, outside existing discovery roots."""
import copy
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import mc
from mc import Material, implicit, explicit, wet_update, yield_value
from column import solve
from compare import META, BULK, CONTACT, compare, load, validate_manifest
from diagnose import tangent_checks, integrate, paths


class LocalChecks(unittest.TestCase):
    def test_closed_form_apex(self):
        m=Material()
        r=implicit(-100*np.eye(3),np.diag([.01,.008,.009]))
        self.assertEqual(r['branch'],'apex')
        expected=m.c/np.tan(np.deg2rad(m.phi))
        np.testing.assert_allclose(r['stress'],expected*np.eye(3),atol=1e-9)

    def test_nonassociated_dilation_and_work(self):
        m=Material()
        de=np.array([[.0003,.003,0],[.003,-.0002,0],[0,0,-.0008]])
        r=implicit(-100*np.eye(3),de,m)
        self.assertEqual(r['branch'],'face')
        lm=float(r['multipliers'].sum())
        self.assertAlmostEqual(np.trace(r['dep']),2*np.sin(np.deg2rad(m.psi))*lm,places=12)
        s=np.linalg.eigvalsh(r['stress'])
        expected=lm*(m.k+(np.sin(np.deg2rad(m.psi))-np.sin(np.deg2rad(m.phi)))*(s[0]+s[-1]))
        self.assertAlmostEqual(r['work'],expected,places=11)
        self.assertGreater(r['work'],0.)

    def test_rotation_covariance(self):
        Q,_=np.linalg.qr(np.random.default_rng(91).normal(size=(3,3)))
        old=np.diag([-60.,-100.,-140.])
        de=np.array([[.0003,.003,0],[.003,-.0002,0],[0,0,-.0008]])
        a=implicit(old,de)
        b=implicit(Q@old@Q.T,Q@de@Q.T)
        np.testing.assert_allclose(b['stress'],Q@a['stress']@Q.T,atol=1e-9)
        np.testing.assert_allclose(b['dep'],Q@a['dep']@Q.T,atol=1e-13)

    def test_return_branches_and_tangent(self):
        rows=tangent_checks()
        self.assertEqual({r['branch'] for r in rows},{'elastic','face','edge','apex'})
        self.assertTrue(all(r['same_branch'] for r in rows))
        # Numerical derivative self-check, NOT an original C-03 acceptance gate.
        self.assertLess(max(r['error_divided_by_E'] for r in rows),1e-7)

    def test_yield_and_nonnegative_work(self):
        rng=np.random.default_rng(1947)
        for _ in range(60):
            e=rng.normal(size=(3,3))*.002
            e=(e+e.T)/2
            r=implicit(-100*np.eye(3),e)
            self.assertLessEqual(yield_value(r['stress']),1e-7)
            self.assertGreaterEqual(r['work'],-1e-10)

    def test_explicit_does_not_call_implicit(self):
        with patch.object(mc,'implicit',side_effect=AssertionError('shared return')):
            r=explicit(-100*np.eye(3),np.diag([.002,.0021,-.01]),32)
        self.assertTrue(np.isfinite(r['stress']).all())

    def test_independent_refinement_difference(self):
        p=paths()['simple_shear']
        errors=[]
        for n in (8,16,32):
            a,b=integrate(p,n,'implicit'),integrate(p,n,'explicit')
            errors.append(np.linalg.norm(a['stress']-b['stress']))
        self.assertTrue(errors[2]<errors[1]<errors[0])
        self.assertTrue(np.isfinite(errors).all())

    def test_wet_effective_stress(self):
        de=np.diag([.002,.0021,-.01])
        dry=implicit(-100*np.eye(3),de)
        for z in (-.25,0.,.5,.75,1.):
            r=wet_update(-100*np.eye(3),de,z)
            np.testing.assert_allclose(r['effective'],dry['stress'],atol=1e-12)
            np.testing.assert_allclose(r['total']+r['pore']*np.eye(3),r['effective'],atol=1e-12)

    def test_column_balance_and_local_energy(self):
        for m,c in ((1,1),(2,1),(1,2)):
            h=solve(m,c)
            self.assertEqual(len(h),24)
            self.assertLess(max(x['free_force_residual'] for x in h),1e-7)
            self.assertLess(max(x['all_node_force_sum'] for x in h),1e-10)
            self.assertLess(max(abs(x['local_energy_identity_error']) for x in h),1e-10)
            self.assertGreaterEqual(min(x['minimum_bulk_dissipation_increment'] for x in h),-1e-10)
            self.assertGreaterEqual(min(x['minimum_friction_dissipation_increment'] for x in h),-1e-10)


class ConsumerChecks(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.m=dict(schema='py-c03-history/1',metadata=dict(META,stress_kind='effective',biot_alpha=1.),
                    schedule=[dict(stage='load',increment=i,time=float(i)) for i in (1,2)],bulk_ids=['B0'],contact_ids=['C0'])
        self.a=dict(schema=self.m['schema'],metadata=copy.deepcopy(self.m['metadata']),records=[])
        for row in self.m['schedule']:
            self.a['records'].append(dict(row,kind='bulk',id='B0',stress=[0.]*6,plastic_strain=[0.]*6,pore=7.5))
            self.a['records'].append(dict(row,kind='contact',id='C0',traction=[0.]*3,plastic_slip=[0.]*3))
        self.write('manifest.json',self.m)
        self.write('a.json',self.a)

    def write(self,name,obj):
        path=self.root/name
        path.write_text(json.dumps(obj),encoding='utf-8')
        return path

    def reject(self,mutate):
        b=copy.deepcopy(self.a)
        mutate(b)
        with self.assertRaises((ValueError,TypeError)):
            load(self.write('bad.json',b),self.m)

    def test_equal_is_still_unaccepted(self):
        r=compare(self.root/'manifest.json',self.root/'a.json',self.root/'a.json')
        self.assertIsNone(r['accepted'])
        self.assertTrue(all(v['max_abs']==0. for v in r['fields'].values()))

    def test_reject_unit_sign_order_stress_kind(self):
        for key,value in dict(stress_unit='Pa',stress_sign='compression_positive',pore_sign='negative_compression',
                              strain_shear='engineering',stress_components=['zz','yy','xx','xy','yz','zx'],stress_kind='total',biot_alpha=True).items():
            with self.subTest(key=key):
                self.reject(lambda b:b['metadata'].__setitem__(key,value))

    def test_reject_missing_gp_increment(self):
        self.reject(lambda b:b['records'].pop())
        self.reject(lambda b:b['records'][0].__setitem__('id','B1'))
        self.reject(lambda b:b['records'][0].__setitem__('increment',3))

    def test_reject_duplicate_and_interpolation(self):
        self.reject(lambda b:b['records'].append(copy.deepcopy(b['records'][0])))
        self.reject(lambda b:b['records'][0].__setitem__('time',1.00000000001))

    def test_reject_nonfinite_components(self):
        self.reject(lambda b:b['records'][0].__setitem__('stress',[0.]*5))
        self.reject(lambda b:b['records'][0].__setitem__('pore',True))
        self.reject(lambda b:b['records'][0].__setitem__('pore',-7.5))
        self.reject(lambda b:b['records'][0].__setitem__('pore',float('nan')))

    def test_reject_manifest_gaps_and_ambiguous_metadata(self):
        m=copy.deepcopy(self.m)
        m['schedule'][1]['increment']=3
        with self.assertRaises(ValueError):validate_manifest(m)
        m=copy.deepcopy(self.m)
        m['metadata']['stress_sign']='compression_positive'
        with self.assertRaises(ValueError):validate_manifest(m)

    def test_duplicate_json_keys(self):
        p=self.root/'duplicate.json'
        p.write_text('{"schema":1,"schema":2}',encoding='utf-8')
        with self.assertRaises(ValueError):load(p,self.m)

    def test_csv_and_sidecar(self):
        cols=['kind','stage','increment','time','id','pore']+[f'{f}_{i}' for f,size in (BULK|CONTACT).items() if size>1 for i in range(size)]
        p=self.root/'a.csv'
        with p.open('w',newline='',encoding='utf-8') as file:
            w=csv.DictWriter(file,fieldnames=cols)
            w.writeheader()
            for row in self.a['records']:
                out={k:row[k] for k in ('kind','stage','increment','time','id')}
                fields=BULK if row['kind']=='bulk' else CONTACT
                for field,size in fields.items():
                    if size==1:out[field]=row[field]
                    else:
                        for i in range(size):out[f'{field}_{i}']=row[field][i]
                w.writerow(out)
        self.write('a.csv.meta.json',dict(schema=self.m['schema'],metadata=self.m['metadata']))
        self.assertEqual(load(p,self.m),load(self.root/'a.json',self.m))
        (self.root/'a.csv.meta.json').unlink()
        with self.assertRaises(FileNotFoundError):load(p,self.m)


if __name__=='__main__':
    unittest.main()
