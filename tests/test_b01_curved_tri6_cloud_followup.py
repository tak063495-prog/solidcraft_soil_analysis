"""Bounded contracts and sealed evidence; no new equilibrium calculation in CI."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2027-01-01-b01-curved-tri6-cloud-followup'
_path=sys.path[:]
try:
    sys.path.insert(0,str(HERE))
    from b01follow_common import EDGES, dot, metric, mm, read, rotate
    from b01follow_evidence import verify
    from b01follow_plan import FIELDS, SHAPES, patch, rejection_controls, rotations, specs, target_fields
    from b01follow_precision import decimal_solve
    from b01follow_report import canonical_recovery
finally:
    sys.path[:]=_path


class B01CurvedTri6CloudFollowup(unittest.TestCase):
    def test_fixed_matrix_and_resource_limit(self):
        p=specs()
        self.assertEqual((len(SHAPES),len(rotations()),len(FIELDS),p['new_global_solves']),(8,4,6,192))
        self.assertEqual((p['maximum_elements'] if 'maximum_elements' in p else p['elements_per_patch'],p['dofs_per_patch']),(8,150))
        self.assertFalse(p['b01_complete'])
        with self.assertRaises(ValueError):patch('large')

    def test_distorted_meshes_keep_actual_straight_midpoints(self):
        for shape in SHAPES:
            x,elements,boundary=patch(shape)
            self.assertEqual((len(x),len(elements),len(boundary)),(25,8,16))
            for nodes in elements:
                for e,(a,b) in enumerate(EDGES):
                    np.testing.assert_array_equal(x[nodes[3+e]],(x[nodes[a]]+x[nodes[b]])*.5)
                self.assertGreater(np.linalg.norm(np.cross(x[nodes[1]]-x[nodes[0]],x[nodes[2]]-x[nodes[0]])),0.)

    def test_rotation_maps_work_and_all_six_dofs(self):
        d=np.array([[.2,-.3,.7,.1,.4,-.5]])
        f=np.array([[1.,2.,3.,4.,5.,6.]])
        for name,q in rotations():
            self.assertAlmostEqual(dot(rotate(d,q),rotate(f,q)),dot(d,f),places=13)
            np.testing.assert_allclose(rotate(rotate(d,q),q.T),d,rtol=0,atol=3e-16)
        self.assertEqual(rotations()[-1][0],'general')

    def test_plane_polynomial_is_compatible_with_drilling_and_shear(self):
        x,_,_=patch('interior-0.12')
        values=target_fields(x,np.eye(3))
        membrane=values['membrane'].reshape(-1,6)
        np.testing.assert_array_equal(membrane[:,5],np.full(25,-.00015))
        # theta_x = w_y, theta_y = -w_x for the prescribed quadratic w.
        bending=values['bending'].reshape(-1,6)
        np.testing.assert_allclose(bending[:,3],.0003*x[:,0]-.0007*x[:,1],rtol=0,atol=0)
        np.testing.assert_allclose(bending[:,4],-.001*x[:,0]-.0003*x[:,1],rtol=0,atol=0)

    def test_fsum_and_zero_reference_are_not_relaxed(self):
        self.assertEqual(dot([1e16,1,-1e16],[1,1,1]),1.)
        self.assertEqual(metric([1e-18],[0.])['relative_l2'],1e-18/1e-12)
        with self.assertRaises(ValueError):metric([math.nan],[0.])
        with self.assertRaises(ValueError):metric([1.],[1.,2.])

    def test_recovery_tensors_return_to_material_coordinates(self):
        q=rotations()[-1][1]
        sample=dict(membrane_strain=[1.,2.,6.],curvature_per_m=[1.,2.,6.],
                    membrane_force_kn_m=[1.,2.,3.],bending_moment_kn=[1.,2.,3.],
                    transverse_shear_strain=[4.,5.],transverse_shear_force_kn_m=[4.,5.],drilling_mismatch_rad=.1)
        record=dict(fields=[dict(name='test',recovery=dict(local_axes=q.T.tolist(),samples=[sample],
                    membrane_energy_kn_m=1.,bending_energy_kn_m=2.,transverse_shear_energy_kn_m=3.,drilling_penalty_energy_kn_m=4.))])
        value=canonical_recovery([record],q,'test')
        np.testing.assert_allclose(value['membrane_strain'],[1.,3.,0.,3.,2.,0.,0.,0.,0.],atol=2e-15,rtol=0)
        np.testing.assert_allclose(value['membrane_force_kn_m'],value['membrane_strain'],atol=0,rtol=0)
        np.testing.assert_array_equal(value['energy_components_kn_m'],[1.,2.,3.,4.])

    def test_decimal_reforms_boundary_rhs_without_changing_matrix(self):
        record=dict(free_dofs=[1],boundary_dofs=[0],stiffness=[[2.,-1.],[-1.,2.]],
                    solves=[dict(name='test',displacement_rotation=[2.,0.],applied_force_moment=[0.,1.])])
        before=json.dumps(record,sort_keys=True)
        self.assertEqual([float(v) for v in decimal_solve(record,'test',50)],[2.,1.5])
        self.assertEqual(json.dumps(record,sort_keys=True),before)

    def test_document_suffix_with_real_crlf(self):
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('B-01は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/name;path.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
                self.assertEqual(path.read_text(encoding='utf-8'),text)

    def test_original_records_and_production_path_unchanged(self):
        for path,digest in read(HERE/'prior-record-binding.json').items():
            self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),digest,path)
        probe=(HERE/'b01follow_probe.rs').read_text()
        self.assertIn('soil_fem_3d::plates::{Section, mindlin::Kernel}',probe)
        self.assertNotIn('include!',probe)
        self.assertNotIn('b01follow', (ROOT/'src/plates/mindlin.rs').read_text())

    def test_saved_status_and_rejections_not_promoted(self):
        for name in ('acceptance.json','measurements.json'):
            r=read(HERE/name)
            self.assertIs(r['passed'],False);self.assertIs(r['b01_complete'],False);self.assertIs(r['diagnostic_only'],True)
            self.assertIsNone(r['curved_element_acceptance'])
        r=read(HERE/'measurements.json')
        self.assertEqual((len(r['cases']),r['new_global_solves'],len(r['rejected_controls'])),(32,192,16))
        self.assertTrue(all(not v['accepted'] for v in r['rejected_controls']))
        self.assertEqual(len(rejection_controls()),4)
        self.assertEqual(r['algebraic_zero_counts'],[6])
        self.assertGreater(r['physical_rigid_outside_nullspace_max'],r['physical_rigid_residual_max'])
        self.assertTrue(r['precision_accounting']['all_50_80_digit_roundtrips_identical'])
        self.assertFalse(read(HERE/'unresolved-cause-record.json')['single_cause_fully_identified'])

    def test_sealed_evidence_and_ci_archive_availability(self):
        r=verify(HERE)
        self.assertTrue(r['all_sha256_verified']);self.assertEqual(r['numerical_reruns'],0)
        index=read(HERE/'archive-index.json')
        self.assertIn('b01follow-formal/solve-ledger.json',index['entries'])
        workflow=(ROOT/'.github/workflows/fem-ci.yml').read_text()
        self.assertNotIn('!/interim/2027-01-01-b01-curved-tri6-cloud-followup/',workflow)


if __name__=='__main__':unittest.main()
