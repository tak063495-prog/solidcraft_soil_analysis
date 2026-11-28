"""Reconstruct p14 results from the archived raw observations."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import unittest
from zipfile import ZipFile

import numpy as np
from plans import HERE, P11, P12, P13, ROUTES
from local_lu import F, audit, elimination
from reload_score import run as score_reload

sys.path.insert(0,str(P11))
from alternate_score import CATEGORIES, P10, P7, STAGES, observed, p10_bytes

ROOT = Path(sys.argv.pop(1)) if __name__=="__main__" and len(sys.argv)>1 else HERE


def read(path):
    return json.loads(path.read_text())


class Evidence(unittest.TestCase):
    def test_chain_and_history_transfer(self):
        for sign in ("plus","minus"):
            old = read(P13/"c02"/(sign+".json"))
            folder = ROOT/"c02"
            doc = read(folder/(sign+".json"))
            self.assertEqual(old["consecutive"],doc["start"])
            radius = old["certified_radius_m"]
            previous = None
            for i,row in enumerate(doc["steps"],doc["start"]+1):
                self.assertEqual((row["step"],row["previous_radius_m"]),(i,radius))
                with np.load(folder/f"{sign}-{i}-history.npz") as h:
                    bad = ~h["leaf_stick"]
                    np.testing.assert_array_equal(np.unique(h["parent"][bad]),row["uncertified_stick_children"])
                    self.assertEqual(np.sum(bad),row["uncertain_leaves"])
                    p = h["paired_at_center_kpa"]
                    self.assertTrue((p[0]<=p[1]).all())
                    self.assertTrue((abs(p)<=h["paired_previous_upper_kpa"][:,None]).all())
                    if previous:
                        for child,(nom,enclosed,pair) in previous.items():
                            take = h["parent"]==child
                            np.testing.assert_array_equal(h["past_nominal_kpa"][:,take],nom)
                            np.testing.assert_array_equal(h["past_enclosed_kpa"][:,take],enclosed)
                            np.testing.assert_array_equal(h["paired_previous_upper_kpa"][take],pair)
                    if row["included"]:
                        previous = {k:(h["out_nominal_kpa"][:,h["parent"]==k],
                                       h["out_enclosed_kpa"][:,h["parent"]==k],
                                       h["out_paired_upper_kpa"][h["parent"]==k]) for k in np.unique(h["parent"])}
                bounds = []
                if sum([row["center_residual_kn"][1],row["past_difference_kn"][1]])*(1+1e-12)<row["lower_times_radius_kn"]:
                    bounds.append(sum([row["center_residual_kn"][1],row["past_difference_kn"][1]])/row["lower_kn_m"]*(1+1e-12)**2)
                for key,file in (("inverse_test","residual"),("center_only_inverse_test","center")):
                    inv=row[key]
                    with np.load(folder/f"{sign}-{i}-{file}.npz") as v:
                        eta=np.linalg.norm(abs(v["after_inverse_m"])+v["after_inverse_radius_m"])
                        self.assertLessEqual(eta,inv["center_after_inverse_upper_m"])
                    q=inv["contraction_upper"]
                    image=inv["center_after_inverse_upper_m"]+q*old["radius_m"]
                    self.assertEqual(image,inv["ball_image_upper_m"])
                    self.assertEqual(inv["included"],q<1 and image<old["radius_m"])
                if row["inverse_test"]["included"]:
                    bounds.append(row["inverse_test"]["root_radius_upper_m"])
                self.assertEqual(row["included"],bool(bounds))
                radius = row["certified_radius_m"]
                if bounds:
                    self.assertAlmostEqual(radius,min(bounds),delta=1e-22)
            self.assertFalse(doc["steps"][-1]["included"])
            self.assertEqual(doc["consecutive"],doc["steps"][-1]["step"]-1)
            self.assertEqual(doc["threshold"],.001)
            self.assertIsNone(doc["full_history_relative_upper"])

    def test_two_8_by_8_solves(self):
        orig=read(P13/"nc/direction-check.json")
        report=read(ROOT/"nc/local-lu.json")
        self.assertEqual((report["global_attempts_added"],report["rtol"],report["atol"]),(0,5e-9,1e-14))
        for name,(r,) in orig.items():
            d=report["paths"][name]
            A=[list(map(F,row)) for row in r["matrix"]]
            b=[-F(v) for v in r["rhs_vector"]]
            x,p,perm,terms=elimination(A,b)
            self.assertEqual(d["native"]["dx"][-1],r["native_dx"])
            self.assertEqual(d["permutation"],perm)
            self.assertEqual(d["exact_pivots"],list(map(float,p)))
            self.assertEqual(sum(terms)+p[-1],A[-1][-1])
            for stage in ("native","one_local_refinement"):
                self.assertEqual(audit(A,b,d[stage]["dx"],x)[0],d[stage])
                native_dx=d[stage]["dx"][-1]
                split=d["existing_binary64_schur"]
                self.assertEqual(split[stage]["passes"],abs(split["direction"]-native_dx)<=1e-14+5e-9*abs(native_dx))

    def test_reload_evidence(self):
        folder=ROOT/"c01/reload"
        audit_prefix=read(folder/"prefix-identity.json")
        self.assertEqual(audit_prefix["common_canonical_sha256"],read(P13/"c01/reload/prefix-identity.json")["common_canonical_sha256"])
        self.assertLessEqual(max(audit_prefix["prefix_charges"].values()),1787)
        import tempfile
        with tempfile.TemporaryDirectory() as temp:
            temp=Path(temp)
            (temp/"trace.jsonl.gz").symlink_to((folder/"trace.jsonl.gz").resolve())
            score_reload(temp)
            self.assertEqual(read(temp/"comparison.json"),read(folder/"comparison.json"))

    def test_disjoint_knots_and_charges(self):
        base=read(P12/"c01/alternate/plan.json")
        full=read(P13/"c01/alternate/point26/plan.json")
        selects=[]
        for name in ROUTES:
            folder=ROOT/"c01/alternate"/name
            p=read(folder/"plan.json")
            selects.append(set(p["added"]))
            self.assertEqual(len(selects[-1]),2)
            self.assertEqual(set(p["knots"]),set(base["p11_knots"])|selects[-1])
            for n in (256,512,1024):
                for v in read(folder/str(n)/"charges.json").values():
                    self.assertEqual(v["total"],v["outer"]+v["extra"])
                    self.assertLessEqual(v["total"],1787)
        self.assertFalse(set.intersection(*selects))
        self.assertEqual(set.union(*selects),set(full["added"]))

    def test_reconstructed_22_fields(self):
        for name in ROUTES:
            folder=ROOT/"c01/alternate"/name
            s=read(folder/"score.json")
            passed=0
            for field,gate in s["fields"].items():
                values=[]
                for n in (256,512,1024):
                    sub=folder/str(n)
                    meta=read(sub/"manifest.json")
                    old=read(P10/str(n)/"manifest.json")["fields"][field]
                    with ZipFile(P7/f"c01/first-four-n{n}.zip") as p, ZipFile(P10/str(n)/"prefix-xor.zip") as x, \
                         ZipFile(P10/str(n)/"fields-after-open.zip") as t, ZipFile(sub/"fields-xor-p10.zip") as d:
                        raw=np.frombuffer(p10_bytes(n,old,field,p,x,t),"u1")
                        size=meta["count"]*int(np.prod(old["shape"]))*np.dtype(old["dtype"]).itemsize
                        raw=raw[:size]
                        if field+".bin" in d.namelist():
                            raw=np.bitwise_xor(raw,np.frombuffer(d.read(field+".bin"),"u1"))
                    self.assertEqual(hashlib.sha256(raw).hexdigest(),meta["fields"][field]["sha256"])
                    values.append(raw.view(old["dtype"]).reshape(meta["count"],*old["shape"])[:s["common_coordinates"]])
                a,b,c=values
                if field=="plastic_slip":
                    v=((a-b)*(b-c)).reshape(len(a),32,3).sum(-1)
                    with np.load(folder/"signed-32.npz") as f:
                        np.testing.assert_array_equal(v,f[name])
                ok=((np.array_equal(a,b) and np.array_equal(b,c)) if field in CATEGORIES
                    else observed(*values)["trend"] in ("identical","decreasing_aligned"))
                self.assertEqual(ok,gate["full"]["passes"])
                passed+=ok
            self.assertEqual(len(s["fields"]),22)
            self.assertEqual(passed,s["passed_fields"])
            self.assertEqual(s["passed"],passed==22 and s["common_coordinates"]==1792)

    def test_matched_recontact_and_following_stages(self):
        from alternate import accepted
        folder=ROOT/"c01/alternate"
        m=read(folder/"comparison.json")
        with gzip.open(folder/"common-recontact.jsonl.gz","rt") as stream:
            rows=[json.loads(line) for line in stream]
        self.assertEqual((len(rows),rows[-1]["t"]),(m["common_recontact_ends"],m["matched_endpoint"]))
        for name in ROUTES:
            dest=folder/name
            clocks=[accepted(dest/str(n)/"trace.jsonl.gz")[1] for n in (256,512,1024)]
            for r in rows:
                a,b=-np.diff([c[r["t"]] for c in clocks],axis=0)
                np.testing.assert_array_equal((a*b).sum(1),r[name])
            with np.load(dest/"signed-32.npz") as f:
                for j,stage in enumerate(STAGES):
                    x=f[name][256*j:256*(j+1)]
                    v=m["full_stage_comparison"][stage][name]
                    if len(x)==256:
                        np.testing.assert_array_equal(v["by_point_m2"],x.sum(0))
                        self.assertEqual(v["sum_m2"],x.sum())
                    else:
                        self.assertIsNone(v)
            if all(1. in c for c in clocks):
                end=read(folder/"at-recontact-end/intervals-7-9.json")["routes"][name]
                np.testing.assert_allclose(end["endpoint_m2"],np.array(m["full_stage_comparison"]["recontact"][name]["endpoint_m2"])[[7,9]],rtol=1e-13,atol=1e-30)
            else:
                self.assertNotIn(name,m["full_endpoint_p11_differences"])

    def test_half_step_accounting_and_stage_totals(self):
        folder=ROOT/"c01/alternate"
        s=read(folder/"intervals-7-9.json")
        m=read(folder/"comparison.json")
        grid=np.array(s["common"])
        bins={k:np.zeros((len(grid)-1,2)) for k in s["routes"]}
        with gzip.open(folder/"signed-halves-7-9.jsonl.gz","rt") as stream:
            for line in stream:
                r=json.loads(line)
                l,h=r["span"]
                i=grid.searchsorted(h)-1
                self.assertTrue(grid[i]<=l<h<=grid[i+1])
                v=(np.array(r["weight_m"])*r["delta_m"]).sum(1)
                np.testing.assert_array_equal(v,r["term_m2"])
                bins[r["path"]][i]+=v
        for name,v in bins.items():
            r=s["routes"][name]
            np.testing.assert_array_equal(v,r["terms_m2"])
            np.testing.assert_allclose(r["initial_m2"]+v.sum(0),r["endpoint_m2"],rtol=1e-8,atol=1e-28)
            source=(folder/name if name in ROUTES else P13/"c01/alternate"/name if name.startswith("point") else P12/"c01/alternate")
            with np.load(source/"signed-32.npz") as f:
                for j,stage in enumerate(STAGES):
                    x=f[name][:m["common_coordinates"]][j*256:(j+1)*256]
                    if len(x):
                        item=m["stages"][stage][name]
                        self.assertEqual(x.sum(),item["sum_m2"])
                        np.testing.assert_array_equal(x[-1],item["endpoint_m2"])


if __name__=="__main__":
    unittest.main(verbosity=2)
