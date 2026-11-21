"""Lightweight arithmetic/transfer regressions; heavy solves use sealed evidence."""
import importlib.util
import math
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
STUDY=ROOT/'interim/2026-12-17-c02-spatial-cause-diagnosis'
sys.path.insert(0,str(STUDY))
import c02diag_transfer as transfer
import c02diag_metrics as metrics
from c02diag_manufactured import manufactured,rule
from c02diag_mapping import ContactPlan,natural


class C02SpatialCauseDiagnosis(unittest.TestCase):
    def test_fixed_reduction_order(self):
        self.assertEqual(transfer.dot([1e16,1.,-1e16],[1.,1.,1.]),1.)
        x=np.array([[3.,4.],[0.,12.]])
        self.assertEqual(transfer.norm(x),13.)
        self.assertEqual(transfer.norm(x,np.ones_like(x)),13.)

    def test_invalid_fields_are_rejected(self):
        for value in ([float('nan')],[float('inf')]):
            with self.assertRaises(ValueError):transfer.norm(value)
        with self.assertRaises(ValueError):transfer.norm([1.],[-1.])
        with self.assertRaises(ValueError):metrics.difference(np.zeros(2),np.zeros(3))
        with self.assertRaises(ValueError):transfer.apply([[1.]],[float('nan')])
        with self.assertRaises(ValueError):transfer.transfer([[0.,0.]],[1.],[[.2,.2]],'unknown')

    def test_projection_reproduces_all_quadratics(self):
        rs,w=rule(4);targets,_=rule(5)
        source=np.asarray([transfer.p2(r) for r in rs]);exact=np.asarray([transfer.p2(r) for r in targets])
        for method in ('p2_l2','cubic'):
            operator=transfer.transfer(rs,w,targets,method)
            self.assertLess(np.max(abs(transfer.apply(operator,source)-exact)),2e-12)

    def test_manufactured_interpolation_orders(self):
        result=manufactured()
        for method,low,high in [('nearest',.8,1.2),('p2_l2',2.8,3.2),('cubic',3.8,4.2)]:
            self.assertEqual(len(result[method]['errors']),4)
            for order in result[method]['orders']:
                self.assertGreater(order,low);self.assertLess(order,high)

    def test_nearest_tie_chooses_lowest_index(self):
        operator=transfer.transfer([[0.,0.],[1.,0.]],[1.,1.],[[.5,0.]],'nearest')
        self.assertEqual(transfer.apply(operator,[7.,9.]).tolist(),[7.])

    def test_time_identity_partition_is_exhaustive(self):
        a0=np.array([[1.,2.,3.],[4.,5.,6.]])
        b0=np.zeros_like(a0);a1=a0*2;b1=a0*.25
        labels=np.array([[[False,True,False],[False,True,True]],[[False,True,True],[False,True,False]]]*2)
        result=metrics.time_identity(a0,b0,a1,b1,labels)
        self.assertEqual(result['identity_roundoff_absolute_l2'],0.)
        self.assertEqual(sum(x['point_times'] for x in result['label_partition'].values()),a0.size)
        self.assertAlmostEqual(sum(x['fraction_change_squared'] for x in result['label_partition'].values()),1.,places=14)

    def test_area_and_point_norms_remain_distinct(self):
        values=[np.array([[0.,0.]]),np.array([[2.,1.]]),np.array([[3.,3.]])]
        point=metrics.triple(values);area=metrics.triple(values,np.array([[1.,9.]]))
        self.assertNotEqual(point['order'],area['order'])
        self.assertEqual(metrics.difference([0.],[0.])['relative_l2'],0.)

    def test_common_mapping_preserves_parent_discontinuity(self):
        from types import SimpleNamespace
        rs,w=rule(4);vertices=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        points=[]
        for root in (0,1):
            for r,weight in zip(rs,w):
                points.append(dict(joint_id='test',face_index=root,nodes_a=[0,1,2],position_m=[r[0],r[1],0.],weight_m2=weight))
        system=SimpleNamespace(points=points)
        targets=np.array([[.2,.3,0.],[.2,.3,0.]])
        for method in ('nearest','cubic','p2_l2'):
            plan=ContactPlan({'nodes':vertices},system,0,targets,np.array([0,1]),method)
            result=plan.apply(np.r_[np.ones(16),np.full(16,9.)])
            self.assertLess(max(abs(result-np.array([1.,9.]))),2e-12)
        with self.assertRaises(ValueError):
            ContactPlan({'nodes':vertices},system,0,np.array([[2.,2.,0.]]),np.array([0]),'nearest')
        with self.assertRaises(ValueError):natural(vertices,np.array([[.2,.3,1.]]))

    def test_quadrature_variants_preserve_affine_patch_and_original(self):
        import json
        sys.path.insert(0,str(ROOT/'tools'))
        from joint_reference import System
        from c02diag_discretization import alternate_system
        mesh=json.loads((ROOT/'examples/joint-water-column-tet10-mesh.json').read_text())
        case=json.loads((ROOT/'examples/stages-joint-water-cycle.json').read_text())
        original=System(mesh,case['joint_system'])
        old_c=original.c.copy();old_points=list(original.points)
        for rule_name,interpolation in [('jacobi4','tri6'),('legendre4','tri6'),('jacobi4','corner3')]:
            system=alternate_system(original,rule_name,interpolation)
            self.assertEqual(len(system.points),32)
            self.assertAlmostEqual(math.fsum(system.weight),1.,places=13)
            u=np.zeros_like(np.asarray(mesh['nodes']))
            for point in system.points:
                for node in point['nodes_b']:u[node,2]=mesh['nodes'][node][2]
            jump=np.einsum('qag,g->qa',system.c,u.ravel())
            exact=np.array([point['position_m'][2] for point in system.points])
            self.assertLess(max(abs(jump[:,2]-exact)),1e-12)
        self.assertTrue(np.array_equal(original.c,old_c))
        self.assertEqual(original.points,old_points)

    def test_hidden_opening_times_are_not_endpoint_labels(self):
        from c02diag_labels import catalog,masks,budget
        coarse=[dict(stage='cycle',progress=t) for t in (.5,1.)]
        fine=[dict(stage='cycle',progress=t) for t in (.25,.5,.75,1.)]
        a=catalog(coarse,np.array([[-1.],[-1.]]))
        b=catalog(fine,np.array([[1.],[-1.],[-1.],[-1.]]))
        mask,duration=masks([a,b],('cycle',),2)
        self.assertEqual(mask.tolist(),[[True],[False]])
        self.assertEqual(duration[0,0],.375)
        report,parts=budget(np.array([[1.],[2.]]),mask,('cycle',),2)
        self.assertEqual(parts['opening_label'].tolist(),[[1.],[1.]])
        self.assertEqual(report['identity_residual_l2'],0.)

    def test_saved_closed_path_primitive_identity(self):
        from c02diag_path import primitives
        jump=np.array([[[0.,0.,-1.]],[[1.,0.,-1.]],[[2.,0.,1.]],[[3.,0.,2.]],[[4.,0.,-1.]]])
        plastic=np.array([[[0.,0.,0.]],[[1.,0.,0.]],[[1.5,0.,0.]],[[1.5,0.,0.]],[[1.5+1/3,0.,0.]]])
        parts=primitives(jump,np.zeros_like(jump),plastic,np.array([[0.,0.,1.]]),np.array([2.]))
        self.assertLess(np.max(abs(parts['defect'])),1e-15)
        self.assertEqual(parts['fraction'][:4,0].tolist(),[1.,1.,.5,0.])

    def test_sealed_evidence_without_heavy_recalculation(self):
        from c02diag_evidence import verify
        result=verify()
        self.assertTrue(result['archive_verified'])
        self.assertFalse(result['c02_passed'])
        self.assertFalse(result['heavy_recalculation_executed'])

    def test_document_status_handles_actual_crlf(self):
        import tempfile
        from c02diag_evidence import document_status
        original=STUDY/'README.ja.md'
        self.assertTrue(document_status(original))
        with tempfile.TemporaryDirectory() as folder:
            copy=Path(folder)/'README.ja.md'
            copy.write_bytes(original.read_text(encoding='utf-8').replace('\n','\r\n').encode('utf-8'))
            self.assertTrue(document_status(copy))

if __name__=='__main__':unittest.main()
