"""Validate measurements and gates from stored observations, independently of solver execution."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import unittest
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
D=Path(sys.argv.pop(1)).resolve() if len(sys.argv)>1 and not sys.argv[1].startswith("-") else HERE
P7=HERE.parent/"2026-11-26-remaining-p7"
P8=HERE.parent/"2026-11-26-remaining-p8"
P9=HERE.parent/"2026-11-26-remaining-p9"
from freeze import freeze_alt,freeze_reload


def read(path):
    return json.loads(path.read_text())


class Evidence(unittest.TestCase):
    def test_preregistered_loads(self):
        actual=read(D/"c01/alternate/plan.json")
        self.assertEqual(freeze_alt(),actual)
        self.assertEqual(len(actual["spans"]),3)
        self.assertEqual(sorted(set(actual["knots"])),actual["knots"])
        self.assertTrue(all(any(s["left"]<=t<=s["right"] for s in actual["spans"]) for t in actual["knots"]))

    def test_reload_freeze(self):
        d=D/"c01/reload"
        self.assertEqual(read(d/"plan.json"),json.loads(json.dumps(freeze_reload(d/"observation-trace.jsonl.gz"))))
        self.assertEqual(hashlib.sha256((d/"plan.json").read_bytes()).hexdigest(),(d/"plan.sha256").read_text().strip())

    def test_endpoint_and_reload_work(self):
        d=D/"c01/reload"
        self.assertEqual(read(d/"prefix-identity.json")["common_canonical_sha256"],read(P9/"c01/reload/prefix-comparison.json")["common_canonical_sha256"])
        z=read(d/"comparison.json")
        self.assertLess(z["controlled_span"],z["previous_span"])
        for m in ("common","p9_coordinate"):
            c=z[m]["counts"]
            for row in c.values():
                self.assertEqual(row["end"],z[m]["t"])
                self.assertEqual(row["attempts"],row["accepted"]+row["rejected"])
                self.assertEqual(row["attempts"]-c["baseline"]["attempts"],row["net_extra_trials"])
        self.assertLessEqual(max(r["attempts"] for r in z["full_budget"].values()),1787)

    def test_matched_contact(self):
        d=D/"c01/alternate"
        rows=read(d/"matched.json")
        self.assertEqual([r["t"] for r in rows],read(d/"plan.json")["knots"]+[1.])
        for row in rows:
            for p in row["points"]:
                a,b=-np.diff(np.array(p["slip"]),axis=0)
                self.assertEqual(float(a@b),p["dot"])
                np.testing.assert_array_equal(a*b,p["components"])
        for j,span in enumerate(read(d/"plan.json")["spans"]):
            statuses=[r["points"][j]["status"] for r in rows if span["left"]<=r["t"]<=span["right"]]
            self.assertTrue(len(statuses)>=5)

    def test_full_field_streams(self):
        d=D/"c01/alternate"
        score=read(d/"score.json")
        self.assertEqual(len(score["fields"]),22)
        self.assertEqual(score["prior_passed"],17)
        for n in (256,512,1024):
            with zipfile.ZipFile(d/str(n)/"fields-after-open.zip") as tail,zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as prefix,zipfile.ZipFile(d/str(n)/"prefix-xor.zip") as delta:
                for name,meta in read(d/str(n)/"manifest.json")["fields"].items():
                    f=name+".bin"
                    size=3*256*int(np.prod(meta["shape"]))*np.dtype(meta["dtype"]).itemsize
                    with prefix.open(f) as a,tail.open(f) as b:
                        data=a.read(size)
                        if f in delta.namelist():
                            data=(np.frombuffer(data,"u1") ^ np.frombuffer(delta.read(f),"u1")).tobytes()
                        sha=hashlib.sha256(data)
                        for chunk in iter(lambda:b.read(1<<20),b""):
                            sha.update(chunk)
                        self.assertEqual(sha.hexdigest(),meta["sha256"])
            self.assertLessEqual(max(r["total"] for r in read(d/str(n)/"charges.json").values()),1787)
        self.assertEqual(score["passed"],score["full_domain"] and score["passed_fields"]==22)

    def test_bounded_fixed_path(self):
        for s in ("plus","minus"):
            b=np.load(D/f"c02/{s}.npz")
            z=np.load(P8/f"c02/{s}/residual.npz")
            self.assertEqual(b["slip_error"].shape,(448,380))
            self.assertTrue(np.all(np.linalg.norm(z["signed"],axis=2)<=b["slip_error"]))
            self.assertTrue(np.all(np.linalg.norm(b["center_slip"]-z["child"][:,:,:2],axis=2)<=b["slip_radius"]))
            summ=read(D/"c02/summary.json")[s]
            self.assertLess(summ["fixed_path_upper"],summ["old_fixed_path_upper"])
            self.assertEqual(summ["phases"]["4"]["reduced_anchors"],64)
            self.assertGreater(summ["p8_field_relative"],.001)

    def test_history_steps(self):
        for s in ("plus","minus"):
            records=[json.loads(l) for l in (D/f"c02/{s}-intervals.jsonl").read_text().splitlines()]
            clock=np.load(P8/f"c02/{s}/kinematics.npz")["clock"]
            np.testing.assert_array_equal([(r["phase"],r["t"]) for r in records],clock[:,:2])
            for r in records:
                if r["history_matrix"] is not None:
                    self.assertGreater(r["coercivity_lower"],0)
                    self.assertTrue(np.all(np.asarray(r["history_matrix"])>=0))
            self.assertTrue(any(r["history_matrix"] is None for r in records))
            self.assertIsNone(read(D/"c02/summary.json")[s]["continuous_reequilibrated_upper"])

    def test_condensed_samples(self):
        for sign in ("plus","minus"):
            z=read(D/f"c02/{sign}-sampled-tangents.json")
            self.assertEqual((z["free"],z["interior"],z["seam"]),(600,450,150))
            self.assertEqual(len(z["samples"]),448)
            self.assertTrue(all(min(r["smallest_singular"])>0 for r in z["samples"]))

    def test_nc_predictions(self):
        z=read(D/"nc/summary.json")
        rows=[json.loads(l) for l in (D/"nc/global.jsonl").read_text().splitlines()]
        self.assertEqual(sum("control" in r for r in rows),z["new_global_calls"])
        for f in z["parents"]:
            j=next(j for j,r in enumerate(rows) if r.get("frozen_prediction",{}).get("parent")==f["parent"])
            later=[r for r in rows[j+1:] if r.get("parent")==f["parent"]]
            for k in ("predicted","ordinary"):
                trials=[r for r in later if r["control"]==k]
                self.assertEqual(trials[0]["width"],f[k])
                self.assertEqual(len(trials),f[k+"_attempts"])
                self.assertTrue(trials[-1]["accepted"])
                self.assertFalse(any(r["accepted"] for r in trials[:-1]))
        self.assertLess(z["load_end"],z["parents"][0]["load"])
        self.assertLess(z["tangent_end"],0)

    def test_nc_partition(self):
        z=read(D/"nc/substeps.json")
        self.assertEqual(len(z["parents"]),12)
        self.assertLessEqual(z["total_global_including_p9"],2000)
        for r in z["parents"]:
            self.assertEqual((r["left"]["count"],r["right"]["count"]),(1,4))
            self.assertLess(r["left"]["width"],r["cut"])
            self.assertGreater(r["right"]["width"],r["cut"])
            for b in r["same_partition_brackets"]:
                i=str(b["substep"])
                self.assertEqual(b["left"]["count"],b["right"]["count"])
                self.assertLess(b["left"]["denominators"][i]*b["right"]["denominators"][i],0)


if __name__=="__main__":
    unittest.main(verbosity=2)
