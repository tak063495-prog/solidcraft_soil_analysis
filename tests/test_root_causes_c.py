"""Regression cases for the three observation-only comparisons."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "interim/2026-11-23-impl-root-causes-c"


def module(name):
    spec = importlib.util.spec_from_file_location("root_causes_" + name, HERE / f"{name}.py")
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


class TestRootCausesC(unittest.TestCase):
    def test_c01_gap_location_and_rejected_only_ablation(self):
        c01 = module("c01")
        gaps = [np.array([-1., 2.]), np.array([1., 2.]), np.array([-1., 2.])]
        observe = c01.proposal_observer(lambda *_: .125, False)
        audit = dict(actual_increment=.5, reject=True)
        self.assertEqual(observe(audit, gaps, .75, 32, 1e-8), .25)
        self.assertEqual(audit["first_gap_zero"], (.375, 0))
        audit["reject"] = False
        self.assertEqual(observe(audit, gaps, .75, 32, 1e-8), .125)

    def test_c01_common_prefix_counts_overshooting_rejection(self):
        c01 = module("c01")
        rows = [dict(start=0, end=.75, accepted=False, interval="crossing"),
                dict(start=0, end=.5, accepted=True, interval="smooth")]
        runs = [dict(route=r, steps=512, trials=dict(open=list(rows),reload=[])) for r in c01.ROUTES]
        runs[0]["trials"]["open"].append(dict(start=.5,end=1.,accepted=True,interval="smooth"))
        result = c01.compare(runs)["gap_locator/open"]
        self.assertEqual(result["common_end"], .5)
        self.assertEqual(result["routes"]["event4"]["crossing"]["rejected"], 1)
        self.assertEqual(result["routes"]["event4"]["smooth"]["accepted_length_sum"], .5)

    def test_c02_alignment_uses_exact_accepted_time(self):
        c02 = module("c02")
        rows = [dict(stage=s, progress=t) for s in c02.STAGES for t in (.25,.5,1.)]
        keys, indices = c02.common_indices(rows, 2)
        self.assertEqual([rows[i]["progress"] for i in indices], [t for _,t in keys])
        for invalid in (rows[:-1], [*rows,rows[-1]], list(reversed(rows)),
                        [*rows[:-1],dict(stage=rows[-1]["stage"],progress=np.nextafter(1.,0))]):
            with self.assertRaises(ValueError):
                c02.common_indices(invalid, 2)

    def test_c02_signed_surface_cancellation_keeps_local_difference(self):
        c02 = module("c02")
        left = np.array([[2.,-2.]])
        result = c02.pair_metrics(left, np.zeros_like(left), np.array([.5,.5]),
                                  [("open",1.)],np.zeros((2,3)),np.array([5,7]))
        self.assertEqual(result["surface_integral"]["absolute_l2"], 0.)
        self.assertGreater(result["area_weighted"]["absolute_l2"], 0.)
        self.assertEqual(result["maximum_at"]["parent_face"], 5)

    def test_c03_requires_only_selected_axis_to_change(self):
        c03 = module("c03")
        def case(m,c):
            return dict(settings=dict(material_maximum_substep=m),joint_system=dict(maximum_increment_m=c))
        self.assertEqual(c03.verify_axis(case(1,2),case(1,3),"contact")[0]["material"],1)
        for right in (case(2,3),case(2,2),case(1,2)):
            with self.assertRaises(ValueError):
                c03.verify_axis(case(1,2),right,"contact")

    def test_c03_joins_classification_to_correct_support(self):
        c03 = module("c03")
        left = dict(u=np.arange(12.).reshape(1,4,3), traction=np.array([[[3.,2.,1.]]]),
                    stress=np.ones((1,2,6)), contact_status=np.array([[2]]),
                    material_return=np.array([[0,1]]))
        right = {k:np.zeros_like(v) for k,v in left.items()}
        material = [dict(element_index=0,local_index=i,shape_functions=[.25]*4,
                         position_m=[0.,0.,0.]) for i in range(2)]
        contacts = [dict(owner_a=0,owner_b=1,joint_id="joint",face_index=3,position_m=[0.,0.,0.],
                         shape=[.5,.5],nodes_a=[0,1],nodes_b=[2,3])]
        contact, bulk = c03.mismatch_rows([dict(stage="open",progress=.5)],left,right,
                                         dict(tetrahedra=[[0,1,2,3]]),material,contacts)
        self.assertEqual(contact["delta_traction_kpa"],[3.,2.,1.])
        self.assertEqual(contact["delta_u_a_m"],[1.5,2.5,3.5])
        self.assertEqual(bulk["gp"],1)
        self.assertEqual(bulk["delta_u_m"],[4.5,5.5,6.5])
        self.assertEqual(bulk["adjacent_traction"][0]["face"],3)
