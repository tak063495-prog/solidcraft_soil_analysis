import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from nc_collapse_bounds import ellipse_support, ellipse_metric, cone_ray, yield_values, discrete_bound, mechanism_bound
from nc_collapse_classify import matrix_diagnostic


class CollapseBoundsTests(unittest.TestCase):
    def test_ellipse_support_duality_and_surface(self):
        rng = np.random.default_rng(5701)
        for _ in range(20):
            d = rng.normal(size=6); pc=123.; m=1.4
            s,h = ellipse_support(d,pc,m)
            self.assertAlmostEqual(float(yield_values(s,pc,m,'ellipse')),0.,places=12)
            self.assertAlmostEqual(h,float(s@d),places=12)
            center=np.r_[np.full(3,pc/2),np.zeros(3)]
            self.assertTrue(np.allclose(ellipse_metric(m)@(s-center), d*(pc/2)/np.sqrt(d@np.linalg.solve(ellipse_metric(m),d))))

    def test_shear_work_pairing(self):
        s,h=ellipse_support(np.array([0,0,0,1.,0,0]),100.,1.4)
        self.assertAlmostEqual(s[3],100*1.4/(2*np.sqrt(3)))
        self.assertAlmostEqual(h,s[3])

    def test_cone_section_support(self):
        for d in np.vstack((np.eye(6),-np.eye(6))):
            s,h=cone_ray(d,1.4)
            self.assertAlmostEqual(float(s[:3].mean()),1.)
            self.assertLessEqual(float(yield_values(s,100.,1.4,'cone')),1e-12)
            self.assertAlmostEqual(h,float(d@s))

    def test_known_single_point_ellipse_bound(self):
        b=np.zeros((1,6,1));b[0,2,0]=-1
        out=discrete_bound(b,[1.],[-100.],[200.],[100.],[1.4])
        expected=(100+50-50*np.sqrt(1+4*1.4**2/9))/200
        self.assertAlmostEqual(out['mu_feasible'],expected,places=7)
        self.assertLess(out['objective_gap'],2e-7)
        self.assertLess(out['equilibrium_relative'],1e-9)
        self.assertLess(out['maximum_yield'],1e-9)
        self.assertFalse(out['physical_bound_certified'])

    def test_axisymmetric_cone_known_formula(self):
        b=-np.eye(6)[None,:,:]
        out=discrete_bound(b,[1.],[100.,100.,0.,0.,0.,0.],[0.,0.,100.,0.,0.,0.],[100.],[1.4],kind='cone')
        self.assertAlmostEqual(out['mu_feasible'],(3-2*1.4)/(3+1.4),places=6)
        self.assertLess(out['equilibrium_relative'],1e-8)
        self.assertLess(out['maximum_yield'],1e-7)

    def test_work_bound_sign_and_witness(self):
        b=np.zeros((1,6,1));b[0,2,0]=-1
        out=mechanism_bound(b,[1.],[-1.],[-100.],[200.],[100.],[1.4])
        self.assertIsNotNone(out['lambda_UB'])
        self.assertEqual(out['admissible_inequality'],'mu >= threshold')

    def test_invalid_data_rejected(self):
        with self.assertRaises(ValueError):ellipse_support(np.ones(6),-1,1.4)
        with self.assertRaises(ValueError):ellipse_support(np.ones(6),1,0)
        with self.assertRaises(ValueError):discrete_bound(np.ones((1,6,1)),[-1],[0],[1],[1],[1])


class MatrixTests(unittest.TestCase):
    def test_regular_bordered_limit_point(self):
        out=matrix_diagnostic(np.diag([1.,0.]),[0.,1.],[0.,1.,0.])
        self.assertEqual(out['k_numerical_rank'],1)
        self.assertEqual(out['augmented_rank'],3)
        self.assertEqual(abs(out['left_load_cosine']),1.)
        self.assertFalse(out['augmented_exact_right_null_exists'])
        self.assertLess(out['rectangular_null_relative_residual'],1e-12)

    def test_bifurcation_orthogonal_load(self):
        out=matrix_diagnostic(np.diag([1.,0.]),[1.,0.],[1.,0.,1.])
        self.assertEqual(out['left_load_cosine'],0.)
        self.assertEqual(out['rectangular_rank'],1)
        self.assertTrue(out['augmented_exact_right_null_exists'])

    def test_nonsymmetric_left_controls_solvability(self):
        out=matrix_diagnostic(np.array([[1.,1.],[0.,0.]]),[1.,0.],[0.,1.,0.])
        self.assertAlmostEqual(out['left_load_cosine'],0.)
        self.assertGreater(abs(out['right_load_cosine']),.5)




class SweepMaterialTests(unittest.TestCase):
    def test_sweep_only_requested_physics_and_separate_cap_pressure(self):
        from nc_collapse_sweep import profile
        from probe_nc_midpoint import inputs
        from copy import deepcopy
        original=inputs('midpoint'); saved=deepcopy(original)
        mesh,ground,case=profile(1.3,50,True)
        self.assertEqual(case,original[2])
        self.assertEqual(ground['initialization']['geostatic']['seed'],original[1]['initialization']['geostatic']['seed'])
        self.assertEqual(ground['initialization']['geostatic']['base_loads']['surface_loads'][0]['pressure_kpa'],100.)
        self.assertEqual(ground['initialization']['geostatic']['base_loads']['surface_loads'][1]['pressure_kpa'],50.)
        self.assertEqual(sorted(mesh['face_sets']['exterior']),sorted(mesh['face_sets']['collapse_caps']+mesh['face_sets']['collapse_sides']))
        self.assertEqual(original,saved)

    def test_gravity_off_is_literal_zero_without_small_weight_substitution(self):
        from nc_collapse_sweep import profile
        _,ground,_=profile(1.4,100,False)
        self.assertTrue(all(m['unit_weight_kn_m3']==0. for m in ground['materials']))
        with self.assertRaises(ValueError):profile(float('nan'),100,True)
        with self.assertRaises(ValueError):profile(1.4,75,True)

    def test_material_request_retains_history_and_engineering_strain(self):
        from nc_collapse_material import request, run
        state={'stress':[1.,2.,3.]};par={'integration':'midpoint'}
        r=request(7,par,state,np.arange(6)/100000.,'test')
        self.assertEqual(r['state'],state);self.assertEqual(r['partitions'],1)
        self.assertEqual(r['strain'],(np.arange(6)/100000.).tolist())
        with self.assertRaises(ValueError):run(None,None,None,strain_steps=0)


class CollapseAppendixTests(unittest.TestCase):
    def setUp(self):
        import importlib.util
        import json
        import tempfile
        root=Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('collapse_appendix_audit',root/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py')
        self.audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.audit)
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.old=json.loads((root/self.audit.PHYSICAL_APPENDICES_PATH).read_text())['appendices']
        self.new=json.loads((root/self.audit.COLLAPSE_APPENDICES_PATH).read_text())['appendices']
        for name in (self.audit.PHYSICAL_APPENDICES_PATH,self.audit.COLLAPSE_APPENDICES_PATH):
            p=self.root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((root/name).read_bytes())
        for name in self.old:
            # The audit pins UTF-8 bytes, including LF. Windows text output
            # translates LF to CRLF even when Git checkout preserves bytes.
            p=self.root/name;p.parent.mkdir(parents=True,exist_ok=True)
            p.write_bytes(('original prefix\n'+self.old[name]+self.new[name]).encode('utf-8'))

    def test_both_pinned_suffixes_checked_before_old_prefix(self):
        self.assertEqual(self.audit.physical_document_appendices(self.root),{p:self.old[p]+self.new[p] for p in self.old})

    def test_modified_new_suffix_is_rejected(self):
        p=self.root/'docs/LOW_PRELOAD_NC.ja.md';p.write_bytes(p.read_bytes()+b'forged')
        with self.assertRaisesRegex(ValueError,'suffix differs'):self.audit.physical_document_appendices(self.root)

    def test_resealed_new_manifest_is_rejected(self):
        p=self.root/self.audit.COLLAPSE_APPENDICES_PATH;p.write_bytes(p.read_bytes()+b' ')
        with self.assertRaisesRegex(ValueError,'SHA-256 differs'):self.audit.physical_document_appendices(self.root)

    def test_missing_new_manifest_is_rejected(self):
        (self.root/self.audit.COLLAPSE_APPENDICES_PATH).unlink()
        with self.assertRaises(FileNotFoundError):self.audit.physical_document_appendices(self.root)

    def test_old_suffix_mutation_is_still_rejected(self):
        p=self.root/'docs/LOW_PRELOAD_NC.ja.md'
        p.write_bytes(('original prefix\nchanged old suffix'+self.new['docs/LOW_PRELOAD_NC.ja.md']).encode('utf-8'))
        with self.assertRaisesRegex(ValueError,'suffix differs'):self.audit.physical_document_appendices(self.root)

    def test_windows_text_translation_is_rejected_without_normalization(self):
        path='docs/LOW_PRELOAD_NC.ja.md';p=self.root/path
        original=p.read_bytes()
        # Explicit Windows text-mode behavior is reproducible on every OS.
        p.write_text(original.decode('utf-8'),encoding='utf-8',newline='\r\n')
        converted=p.read_bytes()
        self.assertNotEqual(converted,original)
        self.assertEqual(converted.replace(b'\r\n',b'\n'),original)
        with self.assertRaisesRegex(ValueError,'NC collapse document appendix suffix differs'):
            self.audit.physical_document_appendices(self.root)

    def test_old_suffix_crlf_is_rejected_after_exact_new_suffix(self):
        path='docs/LOW_PRELOAD_NC.ja.md';p=self.root/path
        p.write_bytes(b'original prefix\n'+self.old[path].encode('utf-8').replace(b'\n',b'\r\n')
                      +self.new[path].encode('utf-8'))
        self.assertEqual(self.audit.collapse_document_appendices(self.root),self.new)
        with self.assertRaisesRegex(ValueError,'NC physical document appendix suffix differs'):
            self.audit.physical_document_appendices(self.root)


class IntersectionTests(unittest.TestCase):
    def test_intersection_support_and_ordering(self):
        from nc_collapse_bounds import intersection_support
        rng=np.random.default_rng(210)
        for _ in range(30):
            d=rng.normal(size=6);s,h=intersection_support(d,100.,1.4)
            self.assertLessEqual(float(yield_values(s,100.,1.4,'ellipse')),1e-12)
            self.assertLessEqual(float(yield_values(s,100.,1.4,'cone')),1e-11)
            self.assertLessEqual(h,ellipse_support(d,100.,1.4)[1]+1e-10)
        b=-np.eye(6)[None,:,:]
        out=discrete_bound(b,[1.],[50.,50.,0.,0.,0.,0.],[0.,0.,100.,0.,0.,0.],[100.],[1.4],kind='intersection')
        self.assertLess(out['maximum_yield'],1e-8)
        self.assertLess(out['objective_gap'],2e-7)


if __name__ == '__main__':unittest.main()
