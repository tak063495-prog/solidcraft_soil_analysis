"""Check p19 against prior immutable equilibria, trial ledgers and physical gates."""
import csv
import gzip
import json
import os
from pathlib import Path
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
P18 = HERE.with_name("2026-11-28-remaining-p18")
OUT = Path(os.environ.get("P19_RESULTS",HERE))


def read(p):
    return json.loads(p.read_bytes())


def jsonl(p):
    return [json.loads(s) for s in p.read_text().splitlines()]


class Evidence(unittest.TestCase):
    def test_p18_face_replay(self):
        for s,d in (("plus",P18/"c02"),("minus",P18/"c02/minus-328")):
            self.assertEqual(read(OUT/f"c02/baseline/{s}.json"),read(d/f"{s}.json"))

    def test_tight_radius(self):
        base,new = [read(OUT/f"c02/{d}/minus.json") for d in ("baseline","tight")]
        for run in (base,new):
            self.assertEqual((run["consecutive"],run["total"],run["threshold"]),(832,4156,.001))
            self.assertIsNone(run["full_history_relative_upper"])
            a,b = run["steps"][-2:]
            self.assertEqual((a["step"],b["step"]),(832,833))
            self.assertLess(a["inverse_test"]["ball_image_upper_m"],1e-8)
            self.assertGreater(b["inverse_test"]["ball_image_upper_m"],1e-8)
            self.assertEqual(a["certified_radius_m"],b["previous_radius_m"])
        self.assertLess(new["steps"][-2]["certified_radius_m"],base["steps"][-2]["certified_radius_m"])
        for name in ("center_only_inverse_test","integration_tiles"):
            self.assertEqual(base["steps"][-1][name],new["steps"][-1][name])
        with np.load(OUT/"c02/tight/minus-832-leaves.npz") as leaves:
            self.assertEqual(len(leaves["selected_indices"]),len(leaves["block_leaf"]))

    def test_leaf_support_and_excess(self):
        all_groups = read(OUT/"c02/leaf-groups.json")
        for name,v in all_groups.items():
            d,s = name.split("-")
            folder = OUT/"c02"/d
            with gzip.open(folder/f"{s}-833-by-leaf.csv.gz","rt") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len({(r["block"],r["leaf"]) for r in rows}),len(rows))
            with np.load(folder/f"{s}-833-pair.npz") as pair:
                for field,term in (("norm_support_nm","remainder_m"),("direction_support_nm","center_direction_support_m")):
                    np.testing.assert_allclose(sum(float(r[field]) for r in rows),pair[term]*1e9,rtol=1e-12)
            np.testing.assert_allclose(v["image_nm"],v["fixed_nm"]+sum(float(r["total_nm"]) for r in rows),rtol=2e-13)
            self.assertLessEqual(v["image_before_last_rank_nm"],10)
            self.assertGreater(v["image_at_last_rank_nm"],10)

    def test_shadow_band_and_budget(self):
        data = read(OUT/"c01/reload/comparison.json")
        self.assertEqual(data["shadow"]["reasons"]["event"],7)
        self.assertEqual(data["shadow"]["reasons"]["event+slip"],2)
        self.assertEqual(data["shadow"]["reasons"]["slip"],14)
        self.assertEqual(sum(p["actual"] for p in data["shadow"]["points"]),9)
        self.assertEqual(data["retrospective_comparator"],"p13_replay")
        self.assertFalse(data["prior_accepted_landing_table_used"])
        for t,cols in data["common"].items():
            self.assertEqual(len(cols),10)
            for v in cols.values():
                self.assertEqual(v["progress"],float(t))
                self.assertEqual(v["attempts"],v["accepted"]+v["rejected"])
        for stage in data["control"]["stages"].values():
            self.assertLessEqual(stage["charged"],1787)
        self.assertEqual(data["control"]["stages"]["reload"]["charged"],data["budget"]["p19"]["attempts"])

    def test_point19_is_observed(self):
        r = read(OUT/"c01/reload/point19.json")
        self.assertTrue(r["observation_only"])
        self.assertEqual((r["attempt"],r["point"],r["exclusion"]),(877,19,1/65536))
        errs = []
        for forecast in r["trials"]:
            self.assertEqual(max(forecast["endpoints"]),r["start"])
            self.assertTrue(forecast["within_guard"])
            errs.append(abs(forecast["root"]-r["event"]))
        self.assertLess(max(errs[1:]),errs[0])
        self.assertGreater(r["end"]-r["event"],r["exclusion"])

    def test_last_crossings(self):
        with (OUT/"c01/alternate/last-crossings.csv").open() as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len(rows),18)
        gates = read(OUT/"c01/alternate/gates.json")
        self.assertEqual(gates,read(P18/"c01/alternate/gates.json"))
        for r in rows:
            route,k,p = r["route"],r["field"],int(r["point"])
            v = gates[route]
            self.assertEqual((len(v["fields"]),v["passed_fields"],v["common_coordinates"]),(22,17,1792))
            with np.load(P18/f"c01/alternate/{route}-reload.npz") as archived:
                a = archived[k+"_growth"][:,p,1]
                i = np.flatnonzero(a<=0)[-1]
                self.assertEqual(len(a),641)
                self.assertTrue((a[i+1:]>0).all())
                np.testing.assert_array_equal(archived["progress"][i:i+2],[float(r["start"]),float(r["end"])])
                np.testing.assert_allclose(float(r["traction_change"])+float(r["field_change"]),a[i+1]-a[i],rtol=1e-12)

    def test_nc_charges(self):
        rows = jsonl(OUT/"nc/global.jsonl")
        result = read(OUT/"nc/summary.json")
        self.assertEqual((result["rtol"],result["atol"]),(5e-9,1e-14))
        self.assertEqual(result["start"],read(P18/"nc/summary.json")["cumulative"])
        expected = list(range(result["start"]+1,result["cumulative"]+1))
        for event in ("start","verdict","diagnostics"):
            self.assertEqual([r["total"] for r in rows if r["event"]==event],expected)
        for v in [r for r in rows if r["event"]=="verdict"]:
            if v["accepted"]:
                self.assertLessEqual(v["iterations"][-1]["relative_force"],1e-12)
                self.assertLessEqual(abs(v["iterations"][-1]["arc_residual_m"]),1e-12)
            else:
                self.assertEqual(v["end_hash"],v["parent_hash"])
            if v["purpose"]=="reconstruct_named":
                archived = read(P18/"nc"/(v["path"]+"-"+v["width_m"].hex()+".json"))
                self.assertEqual((v["end_hash"],v["iterations"]),(archived["end_hash"],archived["iterations"]))

    def test_nc_directions(self):
        result = read(OUT/"nc/summary.json")
        for d in (OUT/"nc").glob("*-checks.jsonl"):
            for r in jsonl(d):
                self.assertEqual(r["allowance"],1e-14+5e-9*abs(r["refined"][-1]))
                self.assertEqual(r["passes"],abs(r["refined"][-1]-r["schur"])<=r["allowance"])
        for key,r in result["paths"].items():
            with gzip.open(OUT/f"nc/{key}-next.json.gz","rt") as f:
                end = json.load(f)
            self.assertEqual(end["state_hash"],r["end_hash"])
            self.assertEqual(end["load"],r["end_lambda"])
            self.assertLess(r["end_lambda"],r["start_lambda"])
            self.assertFalse(r["reached_one"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
