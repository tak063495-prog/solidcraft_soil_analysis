"""Verify the numerical decisions from the saved raw p20 records."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import unittest
import zipfile

import numpy as np

from original import HERE, P19, STAGES, field_bytes, measure
import sys
sys.path.insert(0,str(P19))
from nc_next import parent_at, prior


OUT=Path(os.environ.get("P20_RESULTS",HERE))


def read(path):
    return json.loads(path.read_bytes())


class Evidence(unittest.TestCase):
    def test_face_prefix_and_stop(self):
        for sign in ("minus","plus"):
            p=read(OUT/f"c02/{sign}.json")
            baseline=read(P19/f"c02/baseline/{sign}.json")
            self.assertEqual([r for r in p["steps"] if r["step"]<=832],baseline["steps"][:-1])
            self.assertEqual((p["total"],p["threshold"]),(4156,.001))
            self.assertEqual(p["consecutive"],p["steps"][-1]["step"]-1)
            self.assertFalse(p["steps"][-1]["included"])
            self.assertIsNone(p["full_history_relative_upper"])
            self.assertGreaterEqual(p["steps"][-1]["inverse_test"]["ball_image_upper_m"],1e-8)

    def test_face_rank_and_contributions(self):
        for sign in ("minus","plus"):
            with gzip.open(P19/f"c02/baseline/{sign}-833-by-leaf.csv.gz","rt") as f:
                rows=list(csv.DictReader(f))
            key=lambda r:(int(r["block"]),int(r["leaf"]))
            selected=[key(r) for r in rows if (r["refined"]=="True")==(sign=="plus")]
            kept=[key(r) for r in rows if r["refined"]=="True"] if sign=="minus" else []
            with np.load(OUT/f"c02/{sign}-833-leaves.npz") as data:
                available=set(map(tuple,data["block_leaf"]))
                self.assertEqual(list(map(tuple,data["block_leaf"][data["selected_indices"]])),
                                 [k for k in kept+selected if k in available])
                self.assertEqual(data["child_tile_remainder_m"].shape[1],64 if sign=="minus" else 256)
            with np.load(OUT/f"c02/{sign}-833-pair-supports.npz") as support, np.load(OUT/f"c02/{sign}-833-pair.npz") as pair:
                old_depth=16 if sign=="minus" else 64
                keys=set(map(tuple,support["block_leaf"]))
                self.assertEqual(int(pair["tiles"]),old_depth*(len(keys)+3*len(keys.intersection(selected))))
            s=read(OUT/"c02/summary.json")[sign]
            self.assertAlmostEqual(s["image_nm"],s["fixed_nm"]+sum(g["after"]["total_nm"] for g in s["groups"].values()),places=10)

    def test_original_common_loads(self):
        out=OUT/"c01/original"
        report=read(out/"gates.json")
        self.assertEqual(report["common_coordinates"],1792)
        with zipfile.ZipFile(out/"fields-xor-p7-alternate1024.zip") as z:
            for key,m in read(out/"manifest.json")["fields"].items():
                raw=np.frombuffer(field_bytes(1024,"alternate",key),"u1")
                if key+".bin" in z.namelist():
                    raw=raw^np.frombuffer(z.read(key+".bin"),"u1")
                self.assertEqual(hashlib.sha256(raw).hexdigest(),m["sha256"])
                vals=[np.frombuffer(b,m["dtype"]).reshape(1792,*m["shape"]) for b in
                      [field_bytes(256,"original",key),field_bytes(512,"original",key),raw]]
                self.assertEqual(measure(key,vals),report["fields"][key]["full"])
                for j,s in enumerate(STAGES):
                    self.assertEqual(measure(key,[v[j*256:(j+1)*256] for v in vals]),report["fields"][key]["stages"][s])
        self.assertEqual(len(report["fields"]),22)

    def test_same_environment(self):
        out=OUT/"c01/original"
        score=read(out/"same-environment/gates.json")
        fine=read(out/"manifest.json")
        for key,m in fine["fields"].items():
            vals=[]
            for n in (256,512,1024):
                path=out/f"same-environment/{n}" if n!=1024 else out
                raw=np.frombuffer(field_bytes(n,"original" if n!=1024 else "alternate",key),"u1")
                fname="fields-xor-p7-original.zip" if n!=1024 else "fields-xor-p7-alternate1024.zip"
                with zipfile.ZipFile(path/fname) as z:
                    if key+".bin" in z.namelist():
                        raw=raw^np.frombuffer(z.read(key+".bin"),"u1")
                self.assertEqual(hashlib.sha256(raw).hexdigest(),read(path/"manifest.json")["fields"][key]["sha256"])
                vals.append(raw.view(m["dtype"]).reshape(1792,*m["shape"]))
            self.assertEqual(measure(key,vals),score["fields"][key]["full"])
            for j,stage in enumerate(STAGES):
                self.assertEqual(measure(key,[v[j*256:(j+1)*256] for v in vals]),score["fields"][key]["stages"][stage])

    def test_observer_and_budget(self):
        out=OUT/"c01/original"
        self.assertEqual(read(out/"result.json"),read(P19/"c01/reload/control/result.json"))
        r=read(out/"comparison.json")
        self.assertEqual((len(r["trials"]),r["guard"],max(r["charges"].values())),(101,1/65536,1783))
        with gzip.open(out/"observations.jsonl.gz","rt") as stream:
            events=list(map(json.loads,stream))
        for p,a,t in zip(events[::2],events[1::2],r["trials"],strict=True):
            self.assertLess(p["previous"],p["parent"])
            self.assertEqual((p["parent"],p["to"]),(a["p20_actual"],a["to"]))
            self.assertEqual(p["p20_predict"],t["attempt"])
            for raw,key in ((p,"predicted"),(a,"actual")):
                ratio=np.max(np.abs(raw["error_m"]),axis=1)/raw["budget_m"]
                np.testing.assert_array_equal(ratio,t[key+"_ratio"])
                self.assertEqual(np.flatnonzero(ratio>1).tolist(),t[key])
            self.assertEqual(t["accepted"],not bool(t["actual"]))
        self.assertEqual(sum(not t["accepted"] for t in r["trials"]),15)

    def test_alternate_endpoint_swaps(self):
        rows=read(OUT/"c01/alternate/factors.json")
        self.assertEqual(len(rows),8)
        for r in rows:
            f,m,u=map(np.asarray,(r["projected_field"],r["length_by_level_pair"],r["direction_by_level_pair"]))
            for mask,x in r["swaps_field_magnitude_direction"].items():
                i,j,k=map(int,mask)
                g=np.sum(f[:,i]*u[:,k],axis=-1)*m[:,j]
                self.assertEqual(float(g[1]-g[0]),x)
            values=r["swaps_field_magnitude_direction"]
            self.assertLess(values["000"],0)
            self.assertGreater(values["111"],0)
            np.testing.assert_allclose(sum(r["contributions"].values()),values["111"]-values["000"],rtol=2e-11)
        self.assertEqual(read(OUT/"c01/alternate/gates.json"),read(P19/"c01/alternate/gates.json"))

    def test_global_trials(self):
        s=read(OUT/"nc/summary.json")
        rows=[json.loads(r) for r in (OUT/"nc/global.jsonl").read_text().splitlines()]
        self.assertEqual((s["start"],s["rtol"],s["atol"]),(1406,5e-9,1e-14))
        self.assertLessEqual(s["cumulative"],2000)
        starts,verdicts=rows[::2],rows[1::2]
        self.assertEqual([a["total"] for a in starts],list(range(1407,s["cumulative"]+1)))
        for key,summary in s["paths"].items():
            parent=parent_at(P19/f"nc/{key}-next.json.gz")
            h,lam=prior.digest(parent),parent[1]
            arc=0.
            pos=up=None
            for a,b in zip(starts,verdicts,strict=True):
                self.assertEqual((a["total"],a["path"]),(b["total"],b["path"]))
                if a["path"]!=key:
                    continue
                self.assertEqual((a["parent_hash"],a["lambda_before"],a["arc_m"]),(h,lam,arc))
                if b["accepted"]:
                    last=b["iterations"][-1]
                    self.assertLessEqual(last["relative_force"],1e-12)
                    self.assertLessEqual(abs(last["arc_residual_m"]),1e-12)
                    arc+=a["width_m"]
                    if pos is None and b["equilibrium_tangent"]["lambda_tangent"]>0:
                        pos=b["total"]
                    if up is None and b["lambda_after"]>lam:
                        up=b["total"]
                else:
                    self.assertEqual((b["end_hash"],b["lambda_after"]),(h,lam))
                self.assertEqual(b["arc_m"],arc)
                h,lam=b["end_hash"],b["lambda_after"]
                for d in b["point351_live"]:
                    self.assertEqual(d["allowance"],1e-14+5e-9*abs(d["direction"]))
                    self.assertEqual(d["matched"],abs(d["direction"]-d["schur"])<=d["allowance"])
            for label,value in (("first_positive_tangent",pos),("first_ascending_interval",up)):
                match=summary[label]
                self.assertEqual(match["total"] if match else None,value)
            self.assertEqual((h,lam,arc),(summary["end_hash"],summary["end_lambda"],summary["arc_m"]))
            self.assertEqual(h,prior.digest(parent_at(OUT/f"nc/{key}-end.json.gz")))


if __name__=="__main__":
    unittest.main(verbosity=2)
