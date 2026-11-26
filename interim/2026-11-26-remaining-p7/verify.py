"""Read the frozen evidence and independently re-evaluate budgets, field gates and local directions."""
import gzip
import hashlib
import json
from pathlib import Path
import unittest
import zipfile

import numpy as np

from contact import CATEGORIES, HERE, STAGES, observed


def read(path):
    return json.loads(path.read_text())


def records(path):
    with gzip.open(path,"rt") as stream:
        return list(map(json.loads,stream))


class Evidence(unittest.TestCase):
    def test_shifted_mesh_and_clock(self):
        cost = read(HERE/"c02/cost.json")
        self.assertEqual([cost[k] for k in ("elements","degree","steps","area_points","accepted")],
                         [12,4,512,2660,4156])
        self.assertEqual(cost["attempts_per_stage"],[512,512,512,512,1024,572,512])
        for face in read(HERE/"c02/tiles.json")["faces"]:
            tiles = np.array(face["leaves"])
            areas = np.linalg.det(tiles[:,1:]-tiles[:,:1])
            self.assertTrue((areas>0).all())
            self.assertEqual(sum(areas),1.)
            self.assertTrue(all(m["retained"]<m["count"] for m in face["marks"]))

    def test_surface_norms(self):
        data = read(HERE/"c02/comparison.json")
        native = HERE.parent/"2026-11-25-remaining-p3/c02/p4/principal.npz"
        shifted = HERE/"c02/principal.npz"
        with np.load(native) as a,np.load(shifted) as b:
            for k in ("parent","position","weight"):
                np.testing.assert_array_equal(a[k],b[k])
            for k in ("slip","signed_area_integral"):
                w = np.sqrt(b["weight"])[None,:,None] if k=="slip" else 1.
                rel = np.linalg.norm((a[k]-b[k])*w)/np.linalg.norm(b[k]*w)
                entry = data["pairs"]["native-shifted"][k]
                self.assertAlmostEqual(rel,entry["relative"],places=16)
                self.assertEqual(rel<=.001,entry["within_width"])
        self.assertEqual(hashlib.sha256(native.read_bytes()).hexdigest(),data["sha256"]["native"])
        self.assertEqual(data["width"],.001)

    def test_route_attempts_and_budget(self):
        for route in ("original","alternate"):
            for n in (256,512,1024):
                folder = HERE/f"c01/{route}{n}"
                rows = records(folder/"trace.jsonl.gz")
                result = read(folder/"result.json")
                self.assertEqual(result["complete"],all(x["last_accepted"]==1 for x in result["stages"]))
                self.assertEqual([s["stage"] for s in result["stages"]],list(STAGES))
                for i,s in enumerate(result["stages"]):
                    count = read(folder/"charges.json")[s["stage"]]
                    self.assertEqual(count["extra"],sum(r.get("stage")==i and "charge" in r for r in rows))
                    ordered = [r for r in rows if r.get("stage")==s["stage"]]
                    self.assertEqual(len(ordered),count["outer"])
                    self.assertEqual(count["outer"]+count["extra"],count["total"])
                    self.assertLessEqual(count["total"],1787)
                    pos = 0.
                    for r in ordered:
                        self.assertEqual(r["from"],pos)
                        if r["accepted"]:
                            pos = r["to"]
                            a = r["audit"]
                            self.assertLessEqual(a["max_slip"],a["slip_budget"])
                            self.assertLessEqual(a["max_work"],a["work_budget"])
                        else:
                            self.assertTrue(r["parent_unchanged"])
                    self.assertEqual(pos,s["last_accepted"] or 0.)

    def test_switching_half_rule(self):
        for route in ("original","alternate"):
            for n in (256,512,1024):
                folder = HERE/f"c01/{route}{n}"
                rows = records(folder/"trace.jsonl.gz")
                half = [r for r in rows if "half" in r]
                self.assertTrue(half)
                attempted = {(r["from"],r["to"]):r for r in rows if r.get("stage")=="recontact"}
                for r in half:
                    self.assertEqual(r["slip_budget"],.5*1e-4/1024**2)
                    self.assertEqual(r["work_budget"],.5*.01/n**2)
                    self.assertEqual(r["outside"],r["slip"]>r["slip_budget"] or r["work"]>r["work_budget"])
                    if r["outside"]:
                        self.assertFalse(attempted[r["from"],r["to"]]["accepted"])

    def test_all_field_scores(self):
        for route in ("original","alternate"):
            summary = read(HERE/f"c01/{route}-score.json")
            times = summary["common_coordinates"]
            folders = [HERE/f"c01/{route}{n}" for n in (256,512,1024)]
            meta = [read(d/"manifest.json") for d in folders]
            self.assertEqual(times,min(m["count"] for m in meta))
            passed = times==7*256
            self.assertEqual(len(summary["fields"]),22)
            for name,entry in summary["fields"].items():
                arrays = []
                for n,d,m in zip((256,512,1024),folders,meta):
                    with zipfile.ZipFile(HERE/f"c01/first-four-n{n}.zip") as prefix, zipfile.ZipFile(d/"fields-tail.zip") as tail:
                        raw = prefix.read(name+".bin")+tail.read(name+".bin")
                    f = m["fields"][name]
                    self.assertEqual(hashlib.sha256(raw).hexdigest(),f["sha256"])
                    arrays.append(np.frombuffer(raw,dtype=f["dtype"]).reshape(m["count"],*f["shape"])[:times])
                a,b,c = arrays
                if name in CATEGORIES:
                    counts = [int((a!=b).sum()),int((b!=c).sum())]
                    self.assertEqual(counts,entry["measured"]["mismatch"])
                    ok = counts==[0,0]
                else:
                    measured = observed(a,b,c)
                    self.assertEqual(measured,entry["measured"])
                    ok = measured["trend"] in ("decreasing_aligned","identical")
                passed &= ok
            self.assertEqual(bool(passed),summary["finite_three_level_gate"])

    def test_nc_pair_and_tangent(self):
        a,b = read(HERE/"nc/summary.json")
        self.assertFalse(a["archive"]["accepted"])
        self.assertTrue(b["archive"]["accepted"])
        self.assertEqual(a["point"],b["point"])
        self.assertNotEqual(a["point"],351)
        self.assertEqual(a["archive"]["parent_lambda"],b["archive"]["parent_lambda"])
        self.assertEqual(a["archive"]["ds_m"],2*b["archive"]["ds_m"])
        self.assertEqual(a["archive"]["attempt"]+1,b["archive"]["attempt"])
        for row in (a,b):
            trace = records(HERE/f'nc/substeps-{row["archive"]["attempt"]}.jsonl.gz')
            for r in (r["record"] for r in trace if r["record"].get("iteration")==0):
                J,res = np.array(r["jac_rows"]),np.array(r["r"])
                schur = J[7,7]-J[7,:7]@np.linalg.solve(J[:7,:7],J[:7,7])
                rhs = -res[7]+J[7,:7]@np.linalg.solve(J[:7,:7],res[:7])
                np.testing.assert_allclose(rhs/schur,r["dx"][7],rtol=5e-9,atol=1e-14)
            if not row["archive"]["accepted"]:
                candidates = [r["record"] for r in trace if "alpha" in r["record"]]
                self.assertEqual(len(candidates),30)
                self.assertTrue(all(r["multiplier"]<0 for r in candidates))
        replay = read(HERE/"nc/replay.json")
        self.assertLessEqual(replay["total_including_inherited"],2000)
        self.assertTrue(replay["pair_replayed_exactly"])
        self.assertEqual(replay["unit480_trials"],0)


if __name__=="__main__":
    unittest.main()
