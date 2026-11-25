"""Read-only consistency checks on the retained results and archived solver outputs."""
import gzip
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
P4 = HERE.parent / "2026-11-25-remaining-p4"
P3 = HERE.parent / "2026-11-25-remaining-p3"


class Evidence(unittest.TestCase):
    def test_c02_matched_histories(self):
        root = HERE/"c02"
        summary = json.loads((root/"comparisons.json").read_bytes())
        self.assertEqual(summary["relative_width"], .001)
        folders = {"native":P3/"c02/p4","d4":root/"d4","d8":root/"d8"}
        costs = [json.loads((p/"cost.json").read_bytes()) for p in folders.values()]
        self.assertTrue(all(row["attempts_per_stage"] == costs[0]["attempts_per_stage"] for row in costs))
        fields = {}
        for name, folder in folders.items():
            path = folder/"principal.npz"
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),summary["input_sha256"][name])
            with np.load(path) as f:
                fields[name] = {"slip":f["slip"]*np.sqrt(f["weight"])[None,:,None],
                                "signed_area_integral":f["signed_area_integral"]}
        for key, scores in summary["pairs"].items():
            a,b = key.split("-")
            for name,r in scores.items():
                # Weight each sampled field independently before differencing.
                left,right=fields[a][name],fields[b][name]
                relative=np.linalg.norm(left-right)/np.linalg.norm(right)
                self.assertAlmostEqual(relative,r["relative"],places=14)
                self.assertEqual(relative<=.001,r["within_width"])

    def test_c01_both_routes_and_fixed_parent_errors(self):
        a = json.loads((HERE/"c01/analysis.json").read_bytes())
        for n in (256,512,1024):
            folder = HERE/"c01"/f"n{n}"
            data = []
            with gzip.open(folder/"trials.jsonl.gz","rt") as stream:
                for row in map(json.loads,stream):
                    post=np.subtract(row["after_c"]["slip"],row["after_f"]["slip"])
                    self.assertEqual(abs(post).max(),row["post_max_slip"])
                    if row["accepted"] and (row["end"]*256).is_integer():
                        data.append(row["after_f"]["slip"])
            for path in (P4/"c01/original-fixed-slip.npz",P3/"c01/fixed/slip-common.npz"):
                with np.load(path) as f:
                    np.testing.assert_array_equal(data,f[f"n{n}"][896:1024])
            audit=json.loads((folder/"audit.json").read_bytes())
            self.assertLessEqual(audit["attempts"]+audit["additional_material_and_global_trials"],2065)
        self.assertGreater(a["common_anchors"]["235"]["actual"]["cosine"],0)
        self.assertLess(a["common_anchors"]["236"]["actual"]["cosine"],0)

    def test_nc_same_parent_and_unit_diagnostics(self):
        root=HERE/"nc"
        replay=json.loads((root/"matched-trials.json").read_bytes())["rows"]
        retries=json.loads((root/"matched-retries.json").read_bytes())
        self.assertEqual([r["accepted"] for r in [replay[-1],*retries]],[False,False,False,True])
        self.assertEqual(len({r["parent_lambda"] for r in [replay[-1],*retries]}),1)
        for name in ("line-search","line-search-107","line-search-108"):
            trace=json.loads((root/f"{name}.json").read_bytes())["trace"]
            k=max(i for i,r in enumerate(trace) if "dx" in r)
            self.assertEqual(trace[k]["x"][7],0)
            self.assertLess(trace[k]["dx"][7],0)
            candidates=trace[k+1:]
            self.assertEqual([r["alpha"] for r in candidates],[2.**-i for i in range(30)])
            self.assertTrue(all(r["multiplier"]<0 for r in candidates))
        result=json.loads((root/"unit-comparison.json").read_bytes())
        self.assertEqual(len(result),34)
        for row in result:
            self.assertEqual(row["subnormal"],row["mode"]=="units448")
            self.assertEqual(row["reason"] is None,row["mode"]=="comparison480")
        self.assertEqual(json.loads((root/"attempts.json").read_bytes())["diagnostic_480_global_calls"],0)


if __name__ == "__main__":
    unittest.main()
