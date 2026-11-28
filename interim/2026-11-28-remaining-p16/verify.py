"""Check the measured gates, reconstructed source terms and trial accounting."""
import gzip
import json
import os
from pathlib import Path
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = Path(os.environ.get("P16_RESULTS",HERE))
P15 = HERE.with_name("2026-11-28-remaining-p15")


def read(file):
    return json.loads(file.read_text())


class Evidence(unittest.TestCase):
    def test_paired_prefix(self):
        for sign in ("plus","minus"):
            history = read(OUT/f"c02/{sign}.json")
            previous = read(P15/f"c02/{sign}.json")
            self.assertEqual(history["threshold"],.001)
            self.assertEqual(history["total"],4156)
            self.assertEqual([r for r in history["steps"] if r["step"]<832],previous["steps"][:-1])
            stopped = history["steps"][-1]
            self.assertFalse(stopped["included"])
            self.assertEqual(history["consecutive"],stopped["step"]-1)
            self.assertIsNone(history["full_history_relative_upper"])
            self.assertEqual(next(r for r in history["steps"] if r["step"]==832)["p15_nodal_past"],
                             previous["steps"][-1]["inverse_test"])

    def test_paired_inverse_action(self):
        for sign in ("plus","minus"):
            for r in read(OUT/f"c02/{sign}.json")["steps"]:
                if "p15_nodal_past" not in r:
                    continue
                tag = f"{sign}-{r['step']}"
                with np.load(OUT/f"c02/{tag}-pair.npz") as p,np.load(OUT/f"c02/{tag}-direct.npz") as d:
                    v = r["inverse_test"]
                    nominal = r["center_only_inverse_test"]
                    base = d["center_after_inverse_m"]+p["mid_action_m"]
                    e = min(p["remainder_m"],np.linalg.norm(p["component_radius_m"]))
                    projected = np.sqrt(np.linalg.norm(base)**2 + e**2 +
                        2*(np.linalg.norm(d["center_after_inverse_m"])*p["center_direction_support_m"] +
                           np.linalg.norm(p["mid_action_m"])*e))
                    eta = min(projected,np.linalg.norm(base)+p["remainder_m"],
                              np.linalg.norm(abs(base)+p["component_radius_m"]))*(1+1e-9)+nominal["direct_center_remainder_m"]+v["force_guard_m"]
                    self.assertEqual(eta,v["center_after_inverse_upper_m"])
                    self.assertEqual(eta+1e-8*v["contraction_upper"],v["ball_image_upper_m"])
                    self.assertEqual(v["included"],v["ball_image_upper_m"]<1e-8 and v["contraction_upper"]<1)
                    self.assertGreater(p["tiles"],0)

    def test_original_reject_cost(self):
        b = read(OUT/"c01/reload/p15-bands.json")
        for key in ("width_1_over_16","width_1_over_1024"):
            self.assertEqual(sum(row.get("accepted",0) for row in b[key].values()),1627)
            self.assertEqual(sum(sum(row.values())-row.get("accepted",0) for row in b[key].values()),99)
        top = b["width_1_over_16"]["10"]
        self.assertEqual(top,dict(accepted=343,event=19,slip=11,**{"event+slip":5}))

    def test_proposals_and_comparators(self):
        from reload_cost import run
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            scratch = Path(d)
            (scratch/"trace.jsonl.gz").symlink_to((OUT/"c01/reload/trace.jsonl.gz").resolve())
            run(scratch)
            self.assertEqual(read(scratch/"comparison.json"),read(OUT/"c01/reload/comparison.json"))
        results = read(OUT/"c01/reload/comparison.json")
        self.assertEqual(results["retrospective_comparator"],"p13_replay")
        self.assertFalse(results["prior_accepted_landing_table_used"])
        self.assertEqual(results["budget"]["p16"]["attempts"],1787)
        before = read(P15/"c01/reload/prefix-identity.json")
        after = read(OUT/"c01/reload/prefix-identity.json")
        self.assertEqual(after["common_canonical_sha256"],before["common_canonical_sha256"])

    def test_norm_contributions(self):
        for route in ("outer","middle"):
            for n in (1792,382):
                with np.load(OUT/f"c01/alternate/{route}-{n}.npz") as f:
                    self.assertEqual(len(f["progress"]),n)
                    for field in ("u","jump","reaction"):
                        target = np.diff(f[field+"_squared_norms"],axis=1).ravel()
                        actual = f[field+"_point_growth"].sum(1)+f[field+"_residual_growth"]
                        np.testing.assert_allclose(actual,target,atol=2e-14*f[field+"_squared_norms"].max(),rtol=1e-9)

    def test_common_ends(self):
        with gzip.open(P15/"c01/alternate/common-middle.jsonl.gz","rt") as stream:
            ref = list(map(json.loads,stream))
        for route,key in (("outer","p15_outer"),("middle","p14_middle")):
            with np.load(OUT/f"c01/alternate/{route}-382.npz") as f:
                np.testing.assert_array_equal(f["progress"],[r["t"] for r in ref])
                np.testing.assert_array_equal(f["plastic_slip_dot"],[r["points_m2"][key] for r in ref])

    def test_archived_gate(self):
        from gates import CATEGORIES
        scores = read(OUT/"c01/alternate/gates.json")
        for score in scores.values():
            self.assertEqual(score["common_coordinates"],1792)
            self.assertEqual(len(score["fields"]),22)
            failed = set()
            for k,r in score["fields"].items():
                m = r["full"]
                expected = m["mismatch"]==[0,0] if k in CATEGORIES else m["trend"] in ("identical","decreasing_aligned")
                self.assertEqual(m["passes"],expected)
                if not expected:
                    failed.add(k)
            self.assertEqual(failed,{"u","jump","reaction","plastic_slip","accumulated"})
            self.assertFalse(score["passed"])

    def test_replay_charges(self):
        runs = read(OUT/"c01/alternate/replay-cost.json")
        self.assertEqual(len(runs),6)
        for name,run in runs.items():
            self.assertTrue(run["trace_identical"] or
                            (OUT/"c01/alternate"/(name+"-replay.jsonl.gz")).is_file())
            for s in run["stages"].values():
                self.assertEqual(s["outer"]+s["extra"],s["total"])
                self.assertLessEqual(s["total"],1787)

    def test_parent_decisions(self):
        nc = read(OUT/"nc/parents.json")
        self.assertEqual((nc["rtol"],nc["atol"]),(5e-9,1e-14))
        self.assertEqual(nc["cumulative"]-nc["start_count"],sum(r["additional_global_trials"] for r in nc["paths"].values()))
        for name,r in nc["paths"].items():
            path = OUT/"nc"/(name+"-local-direction.jsonl")
            d, = map(json.loads,path.read_text().splitlines())
            self.assertLessEqual(abs(d["schur"]-d["refined"][-1]),1e-14+5e-9*abs(d["refined"][-1]))
            live = list(map(json.loads,(OUT/"nc"/(name+"-global-direction.jsonl")).read_text().splitlines()))
            self.assertEqual(len(live),r["global_direction_calls"])
            self.assertGreater(len(live),0)
            for v in live:
                self.assertLessEqual(abs(v["schur"]-v["refined"][-1]),1e-14+5e-9*abs(v["refined"][-1]))
            self.assertEqual(r["point"],351)
            self.assertEqual(r["accepted"],r["parent_updated"])
            self.assertEqual(r["accepted"],r["end_hash"]!=r["parent_hash"])
            self.assertEqual(r["additional_global_trials"],1)


if __name__=="__main__":
    unittest.main(verbosity=2)
