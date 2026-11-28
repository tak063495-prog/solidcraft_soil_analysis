"""Evidence checks for the four p11 experiments (no new solver calls)."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import textwrap
import unittest

import numpy as np

HERE=Path(__file__).resolve().parent
P10=HERE.parent/"2026-11-26-remaining-p10"
D=Path(sys.argv.pop(1)).resolve() if len(sys.argv)>1 and not sys.argv[1].startswith("-") else HERE

def read(p):
    return json.loads(p.read_text())

def lines(p):
    return [json.loads(s) for s in p.read_text().splitlines()]

def trace(p,stage):
    with gzip.open(p,"rt") as f:
        return [r for s in f if s.startswith("{") and (r:=json.loads(s)).get("stage")==stage]

class Evidence(unittest.TestCase):
    def test_reload_preregistration(self):
        from crossings import freeze
        p=P10/"c01/reload"
        x=json.loads(json.dumps(freeze(p/"narrow-trace.jsonl.gz",p/"plan.json")))
        self.assertEqual(x,read(D/"c01/reload/plan.json"))
        self.assertEqual((x["attempts"],x["outside_crossed"],x["rejected_crossed"]),(1787,83,56))

    def test_reload_counts_and_proposals(self):
        from reload_score import summarize
        folder=D/"c01/reload"
        self.assertEqual(summarize(folder),read(folder/"comparison.json"))
        self.assertEqual(len({r["progress"] for r in read(folder/"comparison.json")["common"].values()}),1)

    def test_reload_controls_replay(self):
        base=HERE.parent/"2026-11-26-remaining-p9/c01/reload"
        for name,prior in (("ordinary",base/"baseline-trace.jsonl.gz"),("half",P10/"c01/reload/narrow-trace.jsonl.gz")):
            new=trace(D/f"c01/reload/{name}.jsonl.gz","reload")
            self.assertEqual(new,trace(prior,"reload"))
        self.assertEqual(read(D/"c01/reload/prefix-identity.json")["common_canonical_sha256"],
                         read(P10/"c01/reload/prefix-identity.json")["common_canonical_sha256"])

    def test_alternate_half_loads(self):
        p=read(D/"c01/alternate/plan.json")
        self.assertEqual(p["p10_knots"],read(P10/"c01/alternate/plan.json")["knots"])
        self.assertTrue(set(p["p10_knots"])<=set(p["knots"]))
        for point in (17,23):
            r=p["signed_attribution"][str(point)]
            self.assertAlmostEqual(r["net"]+r["initial"]["endpoint_share"],r["end"]["dot"],delta=1e-30)
            for n in (256,512,1024):
                events=[e for e in p["transitions"] if e["point"]==point and e["n"]==n]
                self.assertEqual([e["transition"] for e in events],[["open","slip"]])
                rows=trace(D/f"c01/alternate/{n}/trace.jsonl.gz","recontact")
                shared={x["to"] for x in rows if x["accepted"]}
                if rows and rows[-1]["to"]>=max(p["knots"]) and any(x["accepted"] and x["to"]>=max(p["knots"]) for x in rows):
                    self.assertTrue(set(p["knots"])<=shared)

    def test_joint_stage_budgets(self):
        for n in (256,512,1024):
            p=D/f"c01/alternate/{n}"
            for s,r in read(p/"charges.json").items():
                self.assertEqual(r["outer"]+r["extra"],r["total"])
                self.assertLessEqual(r["total"],1787)
                self.assertLessEqual(r["total"],(2065*12+11)//12)

    def test_alternate_signed_gate(self):
        out=D/"c01/alternate"
        s=read(out/"score.json")
        z=np.load(out/"signed-32.npz")
        self.assertEqual(z["p11"].shape,(s["common_coordinates"],32))
        for j,name in enumerate(("geostatic","shear","open","recontact","reverse","reload","return")):
            r=read(out/"contributions.json").get(name)
            if r:
                v=z["p11"][j*256:(j+1)*256]
                np.testing.assert_array_equal(r["p11"]["by_point_m2"],v.sum(0))
                self.assertEqual(r["p11"]["endpoint_total"],v[-1].sum())
        self.assertEqual(len(s["fields"]),22)
        self.assertEqual(s["passed_fields"],sum(f["full"]["passes"] for f in s["fields"].values()))
        self.assertEqual(s["passed"],s["common_coordinates"]==1792 and s["passed_fields"]==22)

    def test_archived_22_fields(self):
        import zipfile
        from alternate_score import p10_bytes,P7
        for n in (256,512,1024):
            dst=D/f"c01/alternate/{n}"
            old=P10/f"c01/alternate/{n}"
            manifest=read(dst/"manifest.json")
            with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as prefix, zipfile.ZipFile(old/"prefix-xor.zip") as px, \
                    zipfile.ZipFile(old/"fields-after-open.zip") as tail, zipfile.ZipFile(dst/"fields-xor-p10.zip") as dx:
                for name,m in read(old/"manifest.json")["fields"].items():
                    data=p10_bytes(n,m,name,prefix,px,tail)
                    size=manifest["count"]*int(np.prod(m["shape"]))*np.dtype(m["dtype"]).itemsize
                    raw=np.frombuffer(data[:size],"u1")
                    if name+".bin" in dx.namelist():
                        raw=np.bitwise_xor(raw,np.frombuffer(dx.read(name+".bin"),"u1"))
                    self.assertEqual(hashlib.sha256(raw).hexdigest(),manifest["fields"][name]["sha256"])

    def test_face_inventory(self):
        out=D/"c02"
        for sign in ("plus","minus"):
            comp=read(out/"comparison.json")[sign]
            for r in ("0","1e-08"):
                before=lines(P10/f"c02/{sign}-intervals.jsonl")
                paths=[lines(p/f"{sign}-{r}.jsonl") for p in (out,out/"weighted")]
                counts=comp["radii_m"][r]["total"]
                self.assertEqual(counts["p10"],1872)
                self.assertEqual(len(before),4156)
                self.assertEqual(counts["positive"],sum(bool(a["lower"] or b["lower"]) for a,b in zip(*paths,strict=True)))
                self.assertEqual(counts["includes"],sum(bool(a["inclusion"] or b["inclusion"]) for a,b in zip(*paths,strict=True)))
                self.assertEqual(counts["consecutive"],0)
            self.assertEqual(comp["threshold"],.001)
            self.assertIsNone(comp["reequilibrated_continuous_history_relative_upper"])

    def test_nc_prediction_is_p10(self):
        p10=(P10/"nc.py").read_text()
        p11=(HERE/"nc_paths.py").read_text()
        old=p10[p10.index("                flat ="):p10.index("                print(json.dumps(dict(frozen_prediction")]
        new=p11[p11.index("    flat ="):p11.index("    return row")]
        self.assertEqual(textwrap.dedent(old),textwrap.dedent(new))

    def test_nc_independent_paths(self):
        records=lines(D/"nc/global.jsonl")
        report=read(D/"nc/comparison.json")
        for path in ("predicted","ordinary"):
            prev=report["initial_hash"]
            frozen=None
            for r in records:
                if r.get("path")!=path: continue
                if "frozen_prediction" in r:
                    self.assertEqual(prev,r["parent_hash"])
                    frozen=r["frozen_prediction"]["predicted"]
                elif "accepted_parent" in r:
                    prev=r["state_hash"]
                    frozen=None
                else:
                    self.assertEqual(prev,r["hash"])
                    if path=="predicted":
                        self.assertIsNotNone(frozen)
                        self.assertLessEqual(r["width"],frozen)
            m=report["paths"][path]
            self.assertEqual(m["accepted"]+m["rejected"],m["attempts"])
            self.assertEqual(m["arc"],report["common_arc"])
            self.assertLess(m["delta_lambda"],0)
        self.assertEqual(report["total_including_p10"],1182+sum("accepted" in r for r in records))
        self.assertLessEqual(report["total_including_p10"],2000)

if __name__=="__main__":
    unittest.main()
