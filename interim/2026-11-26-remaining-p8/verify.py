"""Recompute P8 checks on the archived inputs and outputs without starting new solves."""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
import zipfile

import numpy as np
import freeze
from face_score import arrays, masks, score, NATIVE
from freeze import HERE,P7


def read(path):
    return json.loads(path.read_bytes())


def rows(path):
    with gzip.open(path,"rt") as f:
        return [json.loads(line) for line in f if line.startswith("{")]


class Evidence(unittest.TestCase):
    def test_frozen_inputs(self):
        with tempfile.TemporaryDirectory() as temp:
            try:
                freeze.HERE=Path(temp)
                freeze.freeze_contact()
                freeze.freeze_nc()
                for name in ("c01/predictor.json","nc/selection.json"):
                    self.assertEqual((freeze.HERE/name).read_bytes(),(HERE/name).read_bytes())
            finally:
                freeze.HERE=HERE

    def test_paired_face_clocks_and_gate(self):
        native=arrays(NATIVE/"principal.npz")
        clock=np.load(NATIVE/"increments.npz")["phase_t_fullforce_linear_yield_iterations"]
        report=read(HERE/"c02/comparison.json")
        for name,sgn in (("plus",1),("minus",-1)):
            folder=HERE/"c02"/name
            cost=read(folder/"cost.json")
            self.assertEqual((cost["degree"],cost["steps"],cost["elements"],cost["area_points"]),(4,512,12,2660))
            self.assertEqual(cost["shift_m"],[sgn/16,sgn/32])
            np.testing.assert_array_equal(np.load(folder/"kinematics.npz")["clock"][:,:2],clock[:,:2])
            ref=arrays(folder/"principal.npz")
            w=np.sqrt(ref["weight"])[None,:,None]
            np.testing.assert_array_equal(native["position"],ref["position"])
            for k,v in (("slip",w),("signed_area_integral",1.)):
                m=score(native[k],ref[k],v)
                gate=report[name]["whole_slip_within" if k=="slip" else "whole_area_within"]
                self.assertEqual(m["relative"]<=.001,gate)

    def test_leaf_defects(self):
        for name in ("plus","minus"):
            folder=HERE/"c02"/name
            records=arrays(folder/"residual.npz")
            field=arrays(folder/"slip.npz")
            parent=np.einsum("tlqc,lq->tlc",field["slip"].reshape(448,380,7,3),field["weight"].reshape(380,7))
            np.testing.assert_array_equal(records["child"]-parent,records["signed"])
            self.assertTrue(records["mixed"].any())
            s=read(HERE/"c02/comparison.json")[name]["residual"]["all"]["mixed_at_anchor"]
            self.assertEqual(s["leaf_anchors"],int(records["mixed"].sum()))
            self.assertEqual(s["max_leaf_residual"],np.linalg.norm(records["signed"][records["mixed"]],axis=-1).max())

    def test_physical_band_integrals(self):
        folder=HERE/"c02"
        a=arrays(folder/"physical-band-integrals.npz")
        r=read(folder/"physical-bands.json")
        mask=masks(a["position"])
        for k,m in mask.items():
            self.assertEqual(int(m.sum()),r["regions"][k]["points"])
            np.testing.assert_allclose(a["area"][m].sum(),r["regions"][k]["area"],atol=2e-15)
            for direction in ("plus","minus"):
                d=score(a["native_"+k],a[direction+"_"+k])
                np.testing.assert_allclose(d["relative"],r[direction][k]["signed_integral"]["relative"],rtol=2e-14)
        for direction in ("native","plus","minus"):
            np.testing.assert_allclose(sum(a[direction+"_"+k] for k in ("west","east","south","north","interior")),
                a[direction+"_all"],atol=1e-16)

    def test_saved_reload_parents(self):
        folder=HERE/"c01"
        trained=read(folder/"predictor.json")["training_rejections"]
        paired=rows(folder/"paired-observations.jsonl.gz")
        counts=read(folder/"reload-comparison.json")["parent_observations"]
        self.assertEqual(len(paired),163)
        self.assertEqual(counts["accepted"],sum(r["accepted"] for r in paired))
        with tarfile.open(folder/"rejected-parents.tar.gz") as archive:
            self.assertEqual(len(archive.getnames()),163)
            for item,r in zip(trained,paired,strict=True):
                self.assertTrue(r["unchanged"])
                raw=archive.extractfile(r["file"]).read()
                self.assertEqual(hashlib.sha256(raw).hexdigest(),r["parent_sha256"])
                self.assertEqual(json.loads(raw)[:3],[item["attempt"],r["from"],r["nominal"]])
                self.assertLessEqual(r["to"],r["nominal"])
                if r["accepted"]:
                    self.assertFalse(r["audit"]["reject"])
                    self.assertFalse(r["switch_reject"])
                    self.assertIsNone(r["error"])

    def test_budgeted_paths(self):
        folder=HERE/"c01"
        frozen=read(folder/"reload-comparison.json")
        for name in ("baseline","predicted"):
            attempts=[r for r in rows(folder/(name+"-trace.jsonl.gz")) if r.get("stage")=="reload"]
            totals=frozen["paths"][name]["reload"]
            self.assertEqual(len(attempts),totals["attempts"])
            self.assertEqual(sum(not r["accepted"] for r in attempts),totals["rejected"])
            self.assertTrue(all(r["parent_unchanged"] for r in attempts if not r["accepted"]))
            self.assertLessEqual(len(attempts),1787)
            for row in attempts:
                if row["accepted"]:
                    a=row["audit"]
                    self.assertLessEqual(a["max_slip"],a["slip_budget"])
                    self.assertLessEqual(a["max_work"],a["work_budget"])
        self.assertTrue(all(x["total"]<=1787 for x in frozen["prefix_charged"].values()))
        self.assertEqual(frozen["paths"]["baseline"]["reload"]["rejected"],163)

    def test_alternate_signed_history(self):
        folder=HERE/"c01/alternate"
        data=arrays(folder/"signed.npz")
        report=read(folder/"stages-points.json")
        for i,n in enumerate(data["levels"]):
            root=P7/f"c01/alternate{n}"
            with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as pre, zipfile.ZipFile(root/"fields-tail.zip") as tail:
                raw=pre.read("plastic_slip.bin")+tail.read("plastic_slip.bin")
            self.assertEqual(hashlib.sha256(raw).hexdigest(),read(root/"manifest.json")["fields"]["plastic_slip"]["sha256"])
            field=np.frombuffer(raw,"<f8").reshape(7,256,32,3)
            np.testing.assert_array_equal(field[4:]-field[3,-1],data["increments"][i])
            np.testing.assert_array_equal(field[3,-1],data["endpoint"][i])
        for label in ("raw","endpoint_removed"):
            x,y=data[label+"_signed_adjacent"]
            for j,k in enumerate(("reverse","reload","return")):
                tab=report[label]["stages"][k]
                self.assertAlmostEqual(tab["cosine"],float(np.sum(x[j]*y[j])/(np.linalg.norm(x[j])*np.linalg.norm(y[j]))))
                expected=(x[j]*y[j]).sum(axis=(0,2))
                np.testing.assert_allclose([r["dot"] for r in sorted(tab["points"],key=lambda r:r["point"])],expected,atol=1e-29)

    def test_stratified_nc(self):
        folder=HERE/"nc"
        frozen=read(folder/"selection.json")["intervals"]
        measured=read(folder/"frequencies.json")
        trace=[json.loads(line) for line in (folder/"global.jsonl").read_text().splitlines()]
        prior=[json.loads(line) for line in (P7/"nc/continuation.jsonl").read_text().splitlines()]
        for a,b in zip(trace,prior,strict=True):
            self.assertEqual((a["accepted"],a["ds"],a["reason"]),(b["accepted"],b["ds_m"],b.get("reason")))
        materials=rows(folder/"materials.jsonl.gz")
        for s,r in zip(frozen,measured["intervals"],strict=True):
            ids={i for g in s["selected"] for i in range(g["first"],g["last"]+1)}
            self.assertEqual(len(ids),r["trials"])
            self.assertEqual(sum(not u["accepted"] for u in trace if u["attempt"] in ids),r["rejected"])
            for i in (351,272):
                selected=[x for x in materials if x["attempt"] in ids and x["input"]["point"]==i]
                starts=[v for x in selected for v in x["observation"]["steps"]]
                self.assertEqual(len(starts),r[f"{i}_newton_starts"])
                self.assertEqual(sum(q["direction"]<0 for q in starts),r[f"{i}_negative_starts"])
                for q in starts:
                    np.testing.assert_allclose(q["rhs"]/q["denominator"],q["direction"],rtol=5e-9,atol=1e-14)
        self.assertEqual(measured["comparison480"],0)
        self.assertLess(measured["load"],1)


if __name__=="__main__":
    unittest.main()
