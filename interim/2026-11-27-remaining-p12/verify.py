"""Validate the p12 evidence from primary records, with no new acceptance criteria."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import unittest
import zipfile

import numpy as np

from alternate_plan import BASE, HERE, LEVELS, halves
from alternate_score12 import STAGES
from alternate_score import P10, P7, p10_bytes, observed, CATEGORIES

OUT=Path(sys.argv.pop(1)) if len(sys.argv)>1 else HERE
P11=BASE.parents[1]


def read(path):
    return json.loads(path.read_text())


def lines(path):
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt") as f:
        return [json.loads(s) for s in f if s.startswith("{")]


class Evidence(unittest.TestCase):
    def test_sequential_inclusion(self):
        for sign in ("plus","minus"):
            rows=lines(OUT/f"c02/{sign}.jsonl")
            r=read(OUT/f"c02/{sign}-summary.json")
            previous=0.
            for i,row in enumerate(rows):
                self.assertEqual(row["step"],i+1)
                self.assertEqual(previous,row["previous_radius_m"])
                self.assertEqual(row["initial_exact_zero"],i==0)
                if i==0:
                    self.assertEqual(row["past_upper_kn"],0.)
                self.assertGreaterEqual(row["total_upper_kn"],row["center_upper_kn"]+row["past_upper_kn"])
                for t in row["tests"]:
                    self.assertEqual(t["included"],t["radius"]*t["lower_kn_m"]>row["total_upper_kn"])
                previous=row["certified_radius_m"]
                if previous is None:
                    self.assertEqual(i,len(rows)-1)
                    self.assertFalse(any(t["included"] for t in row["tests"]))
                else:
                    self.assertTrue(row["tests"][-1]["included"])
                    self.assertLessEqual(previous,row["tests"][-1]["radius"])
            self.assertEqual(r["consecutive"],sum(x["certified_radius_m"] is not None for x in rows))
            self.assertEqual(r["threshold"],.001)
            if r["consecutive"]<4156:
                self.assertIsNone(r["reequilibrated_full_history_relative_upper"])

    def test_entrance_only(self):
        plan=read(OUT/"c01/reload/plan.json")
        old=read(P11/"c01/reload/plan.json")
        self.assertEqual(plan["bands"],old["bands"])
        self.assertEqual(plan["controls"],[r for r in old["controls"] if r["transition"]==["slip","open"]])
        self.assertTrue(all(r["target"]==r["span"][0] for r in plan["controls"]))
        for key in ("common_canonical_sha256","prefix_charges"):
            self.assertEqual(read(OUT/"c01/reload/prefix-identity.json")[key],read(P11/"c01/reload/prefix-identity.json")[key])

    def test_reload_counts(self):
        out=OUT/"c01/reload"
        cmp=read(out/"comparison.json")
        for name,file in (("ordinary",P11/"c01/reload/ordinary.jsonl.gz"),("p11",P11/"c01/reload/cross.jsonl.gz"),("p12",out/"entrance.jsonl.gz")):
            path=[r for r in lines(file) if r.get("stage")=="reload"]
            self.assertLessEqual(len(path),1787)
            self.assertEqual(len(path),cmp["budget"][name]["attempts"])
            c=cmp["common"][name]
            sub=path[:c["attempts"]]
            self.assertTrue(sub[-1]["accepted"])
            self.assertEqual(sub[-1]["to"],cmp["common"]["ordinary"]["progress"])
            self.assertEqual(c["accepted"],sum(r["accepted"] for r in sub))
            self.assertEqual(c["rejected"],sum(not r["accepted"] for r in sub))

    def test_attribution(self):
        plan=read(OUT/"c01/alternate/plan.json")
        terms=lines(OUT/"c01/alternate/signed-halves.jsonl.gz")
        for point in (14,26):
            p=plan["points"][str(point)]
            coeff=p["initial_m2"]+sum(r["term_m2"] for r in terms if r["point"]==point)
            np.testing.assert_allclose(coeff,p["endpoint_m2"],rtol=1e-8,atol=1e-28)
            self.assertEqual(p["interval_term_m2"],min(r[2] for r in p["shared_interval_terms"]))
            np.testing.assert_allclose(sum(r["term_m2"] for r in p["accepted_halves"]),p["interval_term_m2"],rtol=1e-8,atol=1e-28)
            for r in p["open_to_slip"]:
                self.assertEqual(r["labels"],["open","slip"])

    def test_actual_common_loads(self):
        out=OUT/"c01/alternate"
        matches=lines(out/"matched-32.jsonl.gz")
        by_t={r["t"]:r for r in matches if "p12" in r}
        for t in read(out/"plan.json")["knots"]:
            if t<=max(by_t):
                self.assertIn(t,by_t)
        for n in LEVELS:
            charges=read(out/str(n)/"charges.json")
            self.assertTrue(all(x["total"]==x["outer"]+x["extra"]<=1787 for x in charges.values()))
        x=np.load(out/"signed-32.npz")
        report=read(out/"contributions.json")
        for j,stage in enumerate(STAGES):
            for label in ("p11","p12"):
                part=x[label][j*256:(j+1)*256]
                if len(part):
                    self.assertEqual(part.shape[1],32)
                    np.testing.assert_array_equal(report[stage][label]["endpoint"],part[-1])
                    self.assertEqual(report[stage][label]["sum_m2"],part.sum())

    def test_fields_and_gate(self):
        out=OUT/"c01/alternate"
        scores=read(out/"score.json")
        fields={}
        for n in LEVELS:
            sub=out/str(n)
            old=read(P10/str(n)/"manifest.json")["fields"]
            meta=read(sub/"manifest.json")
            with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as prefix, zipfile.ZipFile(P10/str(n)/"prefix-xor.zip") as px, \
                    zipfile.ZipFile(P10/str(n)/"fields-after-open.zip") as tail, zipfile.ZipFile(sub/"fields-xor-p10.zip") as dx:
                for name,m in old.items():
                    raw=np.frombuffer(p10_bytes(n,m,name,prefix,px,tail),"u1")
                    size=meta["count"]*int(np.prod(m["shape"]))*np.dtype(m["dtype"]).itemsize
                    data=raw[:size]
                    if name+".bin" in dx.namelist():
                        data=np.bitwise_xor(data,np.frombuffer(dx.read(name+".bin"),"u1"))
                    self.assertEqual(hashlib.sha256(data).hexdigest(),meta["fields"][name]["sha256"])
                    fields.setdefault(name,[]).append(data.view(m["dtype"]).reshape((meta["count"],*m["shape"]))[:scores["common_coordinates"]])
        self.assertEqual(len(fields),22)
        for name,(a,b,c) in fields.items():
            ok=np.array_equal(a,b) and np.array_equal(b,c) if name in CATEGORIES else observed(a,b,c)["trend"] in ("identical","decreasing_aligned")
            self.assertEqual(ok,scores["fields"][name]["full"]["passes"])
            if name=="plastic_slip":
                a1,b1=-np.diff([a,b,c],axis=0)
                np.testing.assert_array_equal((a1*b1).reshape(len(a1),32,3).sum(-1),np.load(out/"signed-32.npz")["p12"])
        self.assertEqual(scores["passed"],scores["common_coordinates"]==1792 and scores["passed_fields"]==22)

    def test_nc_parent_and_prediction(self):
        out=OUT/"nc"
        records=lines(out/"global.jsonl")
        old=lines(P11/"nc/global.jsonl")
        report=read(out/"comparison.json")
        for path in ("ordinary","predicted"):
            previous=[r["state_hash"] for r in old if r["path"]==path and "accepted_parent" in r][-1]
            start=report["paths"][path]
            self.assertEqual(start["start_hash"],previous)
            accepted=rejected=0
            row=None
            for r in (r for r in records if r["path"]==path and not r.get("replay")):
                if "frozen_prediction" in r:
                    self.assertIsNone(row)
                    self.assertEqual(r["parent_hash"],previous)
                    row=r["frozen_prediction"]
                    upper=min(row["ordinary"],.9*min([*row["partition_boundaries"].values(),
                              *(b["left"]["width"] for b in row["sign_brackets"])]))
                    unique=all(.9*b["left"]["width"]>=upper for b in row["unresolved_brackets"])
                    self.assertEqual(unique,row["unique_width"])
                    self.assertEqual(row["predicted"],upper if unique else None)
                    width=row["ordinary"] if path=="ordinary" else row["predicted"]
                elif "accepted_parent" in r:
                    previous=r["state_hash"]
                    row=None
                else:
                    self.assertIsNotNone(row)
                    self.assertEqual(previous,r["parent_hash"])
                    self.assertLessEqual(r["width"],width)
                    accepted+=int(r["accepted"])
                    rejected+=int(not r["accepted"])
            self.assertEqual((accepted,rejected),(start["accepted"],start["rejected"]))
            if start["stop"] is None:
                self.assertEqual(start["arc"],read(HERE/"design.json")["nc"]["common_arc_m"])
            else:
                self.assertFalse(row["unique_width"])
                self.assertLess(start["arc"], read(HERE/"design.json")["nc"]["common_arc_m"])
            with gzip.open(out/(path+"-end.json.gz"),"rt") as f:
                end=json.load(f)
            data=[end["active_material_states"],end["load"],end["secant"],end["displacement_m"]]
            self.assertEqual(previous,hashlib.sha256(json.dumps(data,separators=(",",":")).encode()).hexdigest())
        self.assertEqual(report["cumulative_including_replay"],1219+report["initial_stop_attempts"]+sum("accepted" in r for r in records))
        self.assertLessEqual(report["cumulative_including_replay"],2000)


if __name__=="__main__":
    unittest.main(verbosity=2)
