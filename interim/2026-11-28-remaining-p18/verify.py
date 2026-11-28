"""Check raw p18 evidence against the fixed physical gates and charged trial records."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = Path(os.environ.get("P18_RESULTS", HERE))
P17 = HERE.with_name("2026-11-28-remaining-p17")


def read(path):
    return json.loads(path.read_text())


class Evidence(unittest.TestCase):
    def test_ranked_nominal_groups(self):
        for sign,step,count in (("minus",832,82),("plus",833,4230)):
            name = f"c02/{sign}-{step}-leaves.npz"
            with np.load(OUT/name) as data, np.load(P17/name) as old:
                tiles = data["tile_remainder_m"]
                np.testing.assert_array_equal(tiles,old["tile_remainder_m"])
                scores = tiles.sum(axis=1)
                sel = data["selected_indices"]
                np.testing.assert_array_equal(sel,np.argsort(-scores,kind="stable")[:count])
                new = data["child_tile_remainder_m"]
                self.assertEqual(new.shape,(count,64))
                self.assertTrue((new>=0).all())
                exp = scores.sum()-scores[sel].sum()+new.sum()
                np.testing.assert_allclose(exp,data["updated_direct_remainder_m"],rtol=3e-13,atol=0)

    def test_full_history_gate(self):
        for sign in ("plus","minus"):
            job = read(OUT/f"c02/{sign}.json")
            self.assertEqual((job["total"],job["threshold"]),(4156,0.001))
            self.assertIsNone(job["full_history_relative_upper"])
            last = job["steps"][-1]
            self.assertEqual(job["consecutive"],last["step"]-1)
            self.assertTrue(all(r["included"] for r in job["steps"][:-1]))
            self.assertFalse(last["included"])
            self.assertGreater(last["inverse_test"]["ball_image_upper_m"],1e-8)
            self.assertIsNone(last["certified_radius_m"])

    def test_enlarged_negative_group(self):
        folder = OUT/"c02/minus-328"
        r = read(folder/"minus.json")
        included,stopped = r["steps"][-2:]
        self.assertEqual((r["consecutive"],included["step"],stopped["step"]),(832,832,833))
        self.assertTrue(included["included"])
        self.assertLess(included["inverse_test"]["ball_image_upper_m"],1e-8)
        self.assertEqual(stopped["previous_radius_m"],included["certified_radius_m"])
        self.assertFalse(stopped["included"])
        for step in (832,833):
            with np.load(folder/f"minus-{step}-leaves.npz") as d:
                np.testing.assert_array_equal(d["selected_indices"],
                    np.argsort(-d["tile_remainder_m"].sum(axis=1),kind="stable")[:328])

    def test_past_refinement(self):
        with np.load(OUT/"c02/plus-833-pair.npz") as new, np.load(P17/"c02/plus-833-pair.npz") as old:
            self.assertEqual(new["tiles"],4*old["tiles"])
            self.assertTrue(np.isfinite(new["remainder_m"]))
            self.assertGreaterEqual(new["remainder_m"],0)

    def test_comparison(self):
        c = read(OUT/"c01/reload/comparison.json")
        self.assertFalse(c["prior_accepted_landing_table_used"])
        self.assertEqual(c["retrospective_comparator"],"p13_replay")
        for k,row in c["common"].items():
            self.assertEqual(len(row),9)
            for value in row.values():
                if value is None:
                    self.assertLess(c["budget"]["p18"]["progress"],float(k))
                    continue
                self.assertEqual(value["progress"],float(k))
                self.assertEqual(value["attempts"],value["accepted"]+value["rejected"])
        self.assertEqual(c["shadow"]["point_forecast"]["event"],22)
        self.assertEqual(c["shadow"]["point_forecast"]["proposals"],367)
        self.assertLessEqual(c["budget"]["p18"]["attempts"],1787)

    def test_reload_return_join(self):
        summary = read(OUT/"c01/alternate/summary.json")
        for route in ("outer","middle"):
            for n in (256,512,1024):
                info = read(OUT/f"c01/alternate/{route}-{n}.json")
                self.assertEqual(len(info["stages"]),7)
                self.assertTrue(all(x["last_accepted"]==1 and x["attempts"]<=1787 for x in info["stages"]))
                self.assertEqual(len(info["trace_sha256"]),64)
            with np.load(OUT/f"c01/alternate/{route}-reload.npz") as a, np.load(P17/f"c01/alternate/{route}-return.npz") as b:
                self.assertEqual(len(a["progress"]),summary["common_endpoints"])
                self.assertTrue((np.diff(a["progress"])>0).all())
                for field in ("u","jump","reaction"):
                    np.testing.assert_array_equal(a[field+"_difference"][:,-1],b[field+"_difference"][:,0])
                    for i in summary["points"]:
                        g = a[field+"_growth"][:,i,1]
                        self.assertLess(g[0],0)
                        self.assertGreater(g[-1],g[-2])
                        self.assertGreater(g[-2],0)
                        np.testing.assert_allclose(g[-1],b[field+"_growth"][0,i,1],rtol=1e-11,atol=1e-35)
                    for n in (256,512,1024):
                        self.assertEqual(a[f"half_{n}"][-1,1],1)
                        self.assertLess(a[f"half_{n}"][-1,0],1)
        self.assertEqual(read(OUT/"c01/alternate/gates.json"),read(P17/"c01/alternate/gates.json"))

    def test_nc_counts_and_signs(self):
        summary = read(OUT/"nc/summary.json")
        rows = [json.loads(s) for s in (OUT/"nc/global.jsonl").read_text().splitlines()]
        self.assertEqual(summary["start"],1368)
        self.assertEqual((summary["rtol"],summary["atol"]),(5e-9,1e-14))
        starts = [r["total"] for r in rows if r["event"]=="start"]
        verdicts = [r["total"] for r in rows if r["event"]=="verdict"]
        self.assertEqual(starts,list(range(1369,summary["cumulative"]+1)))
        self.assertEqual(starts,verdicts)
        for path,info in summary["paths"].items():
            self.assertLessEqual(info["previous_local_success"]["width"],min(info["tested_widths"]))
            self.assertGreaterEqual(info["previous_local_failure"]["width"],max(info["tested_widths"]))
            for k,call in (("predictor",0),("corrected",1)):
                a,b = info[k]
                self.assertLess(a,b)
                for w,sgn in ((a,1),(b,-1)):
                    row = read(OUT/f"nc/{path}-{w.hex()}.json")
                    data = row["calls"][call]
                    self.assertEqual(row["parent_hash"],info["parent_hash"])
                    self.assertEqual(np.sign(data["multiplier"]),sgn)
                    self.assertTrue(data["matched"])
                    self.assertLessEqual(abs(data["multiplier"]-data["schur"]),data["allowance"])
                    self.assertEqual(row["accepted"],row["end_hash"]!=info["parent_hash"])
            with gzip.open(OUT/f"nc/{path}-accepted.json.gz","rt") as f:
                saved=json.load(f)
            self.assertIn(saved["width_m"],info["accepted_widths"])
            checkpoint=[saved["active_material_states"],saved["load"],saved["secant"],saved["displacement_m"]]
            digest=hashlib.sha256(json.dumps(checkpoint,separators=(",",":"),allow_nan=False).encode()).hexdigest()
            self.assertEqual(digest,saved["state_hash"])
            row=read(OUT/f"nc/{path}-{saved['width_m'].hex()}.json")
            self.assertEqual(digest,row["end_hash"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
