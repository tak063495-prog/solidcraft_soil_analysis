"""Checks on persisted evidence; no nonlinear solves."""
import gzip
import json
from pathlib import Path
import unittest

import numpy as np

from freeze_windows import freeze
from score_reload import load_rows

P=Path(__file__).resolve().parent
P8=P.parent/"2026-11-26-remaining-p8"
P7=P.parent/"2026-11-26-remaining-p7"
def read(path):
    return json.loads(path.read_text())


class Evidence(unittest.TestCase):
    def test_preregistered_windows(self):
        self.assertEqual(json.loads(json.dumps(freeze())),read(P/"windows.json"))
        self.assertEqual(len(read(P/"windows.json")["pairs"]),96)

    def test_area_cover(self):
        for sign in ("plus","minus"):
            z=np.load(P/f"c02/{sign}.npz")
            labels=z["class_code"]
            self.assertEqual(labels.shape,(448,380,4))
            self.assertEqual(z["area"].sum(),1.)
            np.testing.assert_array_equal(((labels==3)*z["area"]).sum(axis=2),z["uncertain_area"])
            self.assertTrue(set(np.unique(labels))=={0,1,2,3})

    def test_fixed_history_bounds(self):
        for sign in ("plus","minus"):
            z=np.load(P/f"c02/{sign}.npz")
            old=np.load(P8/f"c02/{sign}/residual.npz")
            rad=z["slip_radius_known"]+z["slip_radius_uncertain"]
            self.assertTrue(np.all(np.linalg.norm(old["child"][:,:,:2]-z["slip_center_integral"],axis=2)<=rad))
            self.assertTrue(np.all(np.linalg.norm(old["signed"],axis=2)<=z["slip_rule_error_upper"]))
            for field in ("traction","slip"):
                upper=z["uncertain_area"]*z[field+"_variation_uncertain"]/2
                self.assertTrue(np.all(z[field+"_radius_uncertain"]<=upper*(1+1e-13)))

    def test_displacement_pair(self):
        z=np.load(P/"c02/equilibrium.npz")
        plus,minus=[np.load(P8/f"c02/{sign}/kinematics.npz") for sign in ("plus","minus")]
        np.testing.assert_array_equal(plus["clock"],minus["clock"])
        np.testing.assert_array_equal(z["u_plus_minus"],plus["u"]-minus["u"])
        self.assertEqual(z["fixed_path_plus_minus_force"].shape,(2,448,250,3))

    def test_reloaded_progress_and_gates(self):
        folder=P/"c01/reload"
        res=read(folder/"comparison.json")
        common=res["matched"]
        for name in ("baseline","predicted","windows"):
            rows=[r for r in load_rows(folder/f"{name}-trace.jsonl.gz") if r.get("stage")=="reload"]
            for r in rows:
                if not r["accepted"]:
                    self.assertTrue(r["parent_unchanged"])
            k=common["paths"][name]["attempts"]
            self.assertTrue(rows[k-1]["accepted"] and rows[k-1]["to"]==common["t"])
            self.assertEqual(common["paths"][name]["rejected"],sum(not r["accepted"] for r in rows[:k]))
            self.assertEqual(len(rows),res["proposal_audit"][name]["verified"])
            self.assertLessEqual(len(rows),1787)
        self.assertLessEqual(max(read(folder/"prefix-comparison.json")["prefix_charges"].values()),1787)

    def test_alternate_endpoint(self):
        z=np.load(P/"c01/alternate/accepted-halves.npz")
        signed=np.load(P8/"c01/alternate/signed.npz")
        for i,n in enumerate((256,512,1024)):
            v=z[f"recontact_{n}_slip"]
            t=z[f"recontact_{n}_progress"]
            self.assertTrue((np.diff(t)>0).all())
            np.testing.assert_array_equal(v[-1],signed["endpoint"][i,[3,5,13]])
            np.testing.assert_array_equal(v[0],z[f"open_{n}_slip"][-1])
            self.assertEqual(len(t),2*read(P7/f"c01/alternate{n}/result.json")["stages"][3]["accepted"]+1)

    def test_signed_onsets(self):
        z=np.load(P/"c01/alternate/accepted-halves.npz")
        s=read(P/"c01/alternate/signs.json")
        for stage,record in s.items():
            for i,p in enumerate((3,5,13)):
                t=sorted(set.intersection(*[set(z[f"{stage}_{n}_progress"]) for n in (256,512,1024)]))
                vals=[z[f"{stage}_{n}_slip"][np.searchsorted(z[f"{stage}_{n}_progress"],t),i] for n in (256,512,1024)]
                products=(vals[0]-vals[1])*(vals[1]-vals[2])
                for dim,negative in (("dot",products.sum(axis=1)<0),("x",products[:,0]<0),("y",products[:,1]<0)):
                    entry=record["points"][str(p)]["exact_common"][dim]["first"]
                    if negative.any(): self.assertEqual(entry["observed"]["t"],t[np.flatnonzero(negative)[0]])
                    else: self.assertIsNone(entry)

    def test_same_load_observations(self):
        folder=P/"c01/alternate"
        report=read(folder/"matched-signs.json")
        plan=read(folder/"matched-plan.json")
        for n in (256,512,1024):
            self.assertLessEqual(max(report["costs"][str(n)].values()),1787)
        for i,point in enumerate((3,5,13)):
            records=report["points"][str(point)]["rows"]
            self.assertEqual([r["t"] for r in records],plan["grids"][i])
            for r in records:
                a,b,c=[np.array(level["slip"]) for level in r["levels"]]
                np.testing.assert_array_equal(a-b,r["coarse_medium"])
                np.testing.assert_array_equal(b-c,r["medium_fine"])
                self.assertEqual(float((a-b)@(b-c)),r["dot"])
                for level in r["levels"]:
                    l,h=level["span"]
                    self.assertTrue(l<r["t"]<=h and level["force"]<=1e-10)

    def test_nc_first_and_simultaneity(self):
        rows=[json.loads(l) for l in (P/"nc/trials.jsonl").read_text().splitlines()]
        path=[r for r in rows if r["kind"]=="path"]
        self.assertEqual([int(r["trial"]) for r in path],list(range(460,638)))
        self.assertEqual([r["trial"] for r in path if r["negative272"]],["636"])
        with gzip.open(P/"nc/point-inputs.jsonl.gz","rt") as f:
            both={}
            for r in map(json.loads,f):
                obs=r["observation"]
                neg={s["step"] for s in obs["steps"] if s["direction"]<0}
                both.setdefault(r["trial"],{}).setdefault(r["evaluation"],{})[obs["point"]]=neg
        for row in rows:
            actual=sorted([e,s] for e,d in both[row["trial"]].items() for s in d[351]&d[272])
            self.assertEqual(actual,row["same_evaluation_and_substep"])

    def test_nc_brackets_and_budget(self):
        summary=read(P/"nc/summary.json")
        rows=[json.loads(l) for l in (P/"nc/trials.jsonl").read_text().splitlines()]
        table={(r["group"],r["width"]):r for r in rows}
        for b in summary["brackets"]:
            for side in ("lower","upper"):
                self.assertEqual(bool(table[b["group"],b[side]][b["field"]]),b[side+"_sign"])
            self.assertLessEqual(b["upper"]-b["lower"],summary["grid_base"]/16/1024*1.000001)
        supplement=read(P/"nc/partition.json")
        self.assertEqual(supplement["total_charged_conservatively"],summary["total_charged"]+supplement["additional_calls"])
        self.assertLessEqual(supplement["total_charged_conservatively"],2000)


if __name__=="__main__":
    unittest.main()
