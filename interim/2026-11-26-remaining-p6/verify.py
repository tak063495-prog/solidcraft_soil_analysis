"""Read-only checks of the committed measurements and archived parents."""
import gzip
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

P = Path(__file__).resolve().parent
P5 = P.parent/"2026-11-26-remaining-p5"
ROOT = P.parents[1]


def data(path):
    return json.loads(path.read_bytes())


def rows(path):
    with gzip.open(path,"rt") as stream:
        return list(map(json.loads,stream))


class Evidence(unittest.TestCase):
    def test_face_norms_and_archives(self):
        s=data(P/"c02/summary.json")
        self.assertEqual(s["width"],.001)
        for filename,digest in s["sha256"].items():
            self.assertEqual(hashlib.sha256((ROOT/filename).read_bytes()).hexdigest(),digest)
        for pair in s["pairs"].values():
            for item in pair.values():
                self.assertEqual(item["relative"],item["numerator"]/item["denominator"])
                self.assertEqual(item["within_width"],item["relative"]<=.001)
                self.assertAlmostEqual(sum(item["stage_squared_share"]),1.)
        for label in ("local1","focused","local2"):
            d=P/"c02"/label
            with np.load(d/"slip.npz") as raw,np.load(d/"principal.npz") as common:
                np.testing.assert_array_equal(np.einsum("tqc,q->tc",raw["slip"],raw["weight"]),
                                              common["signed_area_integral"])

    def test_face_partition_and_clock(self):
        time_grid=data(P.parent/"2026-11-25-remaining-p3/c02/p4/cost.json")["attempts_per_stage"]
        for label in ("local1","focused","local2"):
            d=P/"c02"/label
            cost=data(d/"cost.json")
            self.assertEqual(cost["attempts_per_stage"],time_grid)
            self.assertLessEqual(max(time_grid),2065)
            self.assertEqual((cost["degree"],cost["steps"],cost["elements"]),(4,512,12))
            faces=data(d/"tiles.json")["faces"]
            self.assertEqual(sum(f["points"] for f in faces),cost["area_points"])
            for f in faces:
                self.assertGreater(f["tiles_by_depth"][0],0)
                leaves=np.array(f["leaves"])
                self.assertAlmostEqual(np.abs(np.linalg.det(leaves[:,1:]-leaves[:,:1])).sum(),1.)

    def test_joint_parents_and_cost(self):
        s=data(P/"c01/summary.json")
        for n in (256,512,1024):
            d=P/f"c01/n{n}"
            a=data(d/"audit.json")
            ref=rows(P5/f"c01/n{n}/trials.jsonl.gz")
            cand=rows(d/"candidates.jsonl.gz")
            parents=rows(d/"parents.jsonl.gz")
            self.assertEqual(a["charged"],a["base_attempts"]+a["p5_extra"]+sum(r["extra_solves"] for r in cand))
            self.assertLessEqual(a["charged"],2065)
            self.assertEqual(a["selected"],[r["index"] for r in cand])
            self.assertEqual(a["selected"],[r["index"] for r in parents])
            self.assertTrue(a["original_and_alternate_prefix_equal"])
            switches=set()
            for c,parent in zip(cand,parents):
                r=ref[c["index"]]
                for key in ("status","slip"):
                    self.assertEqual(c["parent"][key],r["parent"][key])
                self.assertEqual(parent["parent"]["stage_progress"],c["start"])
                self.assertEqual(len(c["fine"]["status"]),32)
                self.assertEqual(c["slip_budget"],1e-4/1024**2)
                self.assertEqual(np.max(np.abs(np.array(c["coarse"]["slip"])-c["fine"]["slip"])),r["post_max_slip"])
                if c["accepted"]:
                    for i,(x,y,z) in enumerate(zip(*(c[k]["status"] for k in ("parent","half","fine")))):
                        if (x,y)==("open","slip") or (y,z)==("open","slip"):
                            switches.add(i)
            self.assertEqual(len(switches),32)
            for item in s[str(n)]["regions"].values():
                self.assertEqual(sum(item["reasons"].values()),item["trials"])

    def test_nc_parent_and_scalar_direction(self):
        parent=json.loads(gzip.decompress((P5/"nc/parent.json.gz").read_bytes()))
        self.assertEqual(parent["load"],.5481398948309102)
        summaries=data(P/"nc/summary.json")
        self.assertEqual([r["archive"]["accepted"] for r in summaries],[False,False,False,True])
        for r in summaries:
            a=r["archive"]["attempt"]
            records=[json.loads(l) for l in (P/f"nc/inputs-{a}.jsonl").read_text().splitlines()]
            self.assertEqual(len(records),len(r["evaluations"]))
            for record in records:
                self.assertEqual(record["state"],parent["active_material_states"][351])
        slopes=data(P/"nc/consistency-slope.json")
        self.assertEqual([s["negative_initial_direction"] for s in slopes],[[23],[23],[],[22],[],[],[]])
        for s in slopes:
            self.assertFalse(s["first_yield_splits"])
            for step in s["steps"]:
                np.testing.assert_allclose(step["rhs"]/step["schur"],step["direction"],rtol=5e-9)

    def test_route_resource_limits(self):
        for name in ("local1","focused","local2","nc","c01-n256","c01-n512","c01-n1024"):
            log=data(P/f"checks/{name}.log")
            self.assertEqual(log["exit_code"],0)
            self.assertLess(log["peak_rss_kib"]*1024,40e9)
            self.assertLess(log["seconds"],86400)


if __name__=="__main__":
    unittest.main()
