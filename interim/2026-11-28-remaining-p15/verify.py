"""Cross-check p15 evidence, including failed gates and every trial charge."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
P14=HERE.with_name("2026-11-27-remaining-p14")
OUT=Path(os.environ.get("P15_RESULTS",HERE))
from trial_cost import cost,load_trace


def read(p):
    return json.loads(p.read_text())


class Evidence(unittest.TestCase):
    def test_c02_prefix_and_first_stop(self):
        for sign in ("plus","minus"):
            old=read(P14/f"c02/{sign}.json")
            new=read(OUT/f"c02/{sign}.json")
            self.assertEqual(old["consecutive"],830)
            self.assertEqual(new["threshold"],.001)
            for i,r in enumerate(new["steps"],new["start"]+1):
                self.assertEqual(r["step"],i)
                if i==831:
                    self.assertEqual(r["nodal_full_history"],old["steps"][-1]["inverse_test"])
                if i>new["start"]+1:
                    self.assertEqual(r["previous_radius_m"],new["steps"][i-new["start"]-2]["certified_radius_m"])
            self.assertEqual(new["consecutive"],new["steps"][-1]["step"]-1)
            self.assertIsNone(new["full_history_relative_upper"])

    def test_c02_inverse_and_inclusion(self):
        for sign in ("plus","minus"):
            for r in read(OUT/f"c02/{sign}.json")["steps"]:
                if "integration_tiles" not in r:
                    continue
                with np.load(OUT/f"c02/{sign}-{r['step']}-direct.npz") as v:
                    self.assertEqual(int(v["tiles"]),r["integration_tiles"])
                    self.assertGreater(v["direct_error_m"],0)
                for key in ("center_only_inverse_test","inverse_test"):
                    t=r[key]
                    image=t["center_after_inverse_upper_m"]+t["contraction_upper"]*1e-8
                    self.assertEqual(image,t["ball_image_upper_m"])
                    self.assertEqual(t["included"],t["contraction_upper"]<1 and image<1e-8)
                self.assertEqual(r["included"],r["certified_radius_m"] is not None)

    def test_reload_trial_log_and_fixed_costs(self):
        prior=read(P14/"c01/reload/prefix-identity.json")
        actual=read(OUT/"c01/reload/prefix-identity.json")
        for key in ("common_canonical_sha256","prefix_charges"):
            self.assertEqual(actual[key],prior[key])
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/"trace.jsonl.gz").symlink_to((OUT/"c01/reload/trace.jsonl.gz").resolve())
            subprocess.run([sys.executable,HERE/"reload_score.py",p],check=True)
            self.assertEqual(read(p/"comparison.json"),read(OUT/"c01/reload/comparison.json"))

    def test_outer_pretrial_gap_history(self):
        for n in (256,512,1024):
            trace=load_trace(OUT/f"c01/alternate/outer/{n}/trace.jsonl.gz")
            trials=[r for r in trace if r.get("stage")=="recontact" and "accepted" in r]
            parents=[]
            for r in (r for r in trace if "outer_parent" in r):
                t=r["outer_parent"]
                trial=trials[r["attempt"]-1]
                self.assertEqual((t,r["to"]),(trial["from"],trial["to"]))
                if not parents or parents[-1][0]!=t:
                    parents.append([t,r["gaps"]])
                prior=parents[-2] if len(parents)>1 else None
                self.assertEqual(prior,r["prior"])
                proposal=r["nominal"]
                if prior:
                    for old,gap in zip(prior[1],r["gaps"],strict=True):
                        if old>gap>0:
                            proposal=min(proposal,t+max(gap*(t-prior[0])/(old-gap),r["grid"]/512))
                self.assertEqual(proposal,r["to"])
            self.assertGreater(len(parents),1)

    def test_charges_and_parent_commit(self):
        for n in (256,512,1024):
            d=OUT/f"c01/alternate/outer/{n}"
            totals,bins,_=cost(d/"trace.jsonl.gz")
            charges=read(d/"charges.json")["recontact"]
            self.assertEqual([totals[k] for k in ("outer","extra","total")],
                             [charges[k] for k in ("outer","extra","total")])
            self.assertLessEqual(totals["total"],1787)
            self.assertEqual(sum(r["rejected"] for r in bins.values()),totals["rejected"])

    def test_22_field_hashes_and_gate(self):
        p11=HERE.with_name("2026-11-26-remaining-p11")
        sys.path.insert(0,str(p11))
        from alternate_score import P7,P10,p10_bytes
        root=OUT/"c01/alternate/outer"
        for n in (256,512,1024):
            src=P10/str(n)
            d=root/str(n)
            old=read(src/"manifest.json")["fields"]
            manifest=read(d/"manifest.json")
            new=manifest["fields"]
            with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as prefix, zipfile.ZipFile(src/"prefix-xor.zip") as xor, \
                    zipfile.ZipFile(src/"fields-after-open.zip") as tail, zipfile.ZipFile(d/"fields-xor-p10.zip") as patches:
                for k,m in new.items():
                    base=p10_bytes(n,old[k],k,prefix,xor,tail)
                    size=int(np.prod(m["shape"]))*np.dtype(m["dtype"]).itemsize*manifest["count"]
                    data=np.frombuffer(base[:size],"u1")
                    if k+".bin" in patches.namelist():
                        delta=np.frombuffer(patches.read(k+".bin"),"u1")
                        self.assertEqual(len(delta),size)
                        data=np.bitwise_xor(data,delta)
                    self.assertEqual(hashlib.sha256(data).hexdigest(),m["sha256"])
        s=read(root/"score.json")
        self.assertEqual(len(s["fields"]),22)
        self.assertEqual(s["passed"],s["common_coordinates"]==1792 and s["passed_fields"]==22)

    def test_schur_compensations_and_unchanged_tolerance(self):
        r=read(OUT/"nc/schur.json")
        self.assertEqual((r["rtol"],r["atol"],r["parent_updates"],r["new_global_trials"]),(5e-9,1e-14,0,0))
        self.assertEqual(set(r["paths"]),{"ordinary","predicted"})
        for v in r["paths"].values():
            self.assertEqual(v["point"],351)
            for variant in v["variants"].values():
                self.assertEqual(variant["passes"],abs(variant["direction"]-v["lu"])<=v["allowance"])
            self.assertFalse(v["variants"]["existing"]["passes"])
            self.assertTrue(v["variants"]["both"]["passes"])


if __name__=="__main__":
    unittest.main(verbosity=2)
