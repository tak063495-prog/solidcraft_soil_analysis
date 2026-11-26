"""Check p13 results against committed traces, checkpoints and reconstructed fields."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import unittest
import zipfile

import numpy as np
from plans import HERE, P11, P12
sys.path.insert(0,str(P11))
from alternate_score import P10, P7, STAGES, CATEGORIES, p10_bytes, observed

OUT=Path(sys.argv.pop(1)) if __name__=="__main__" and len(sys.argv)>1 else HERE


def read(path):
    return json.loads(path.read_text())


def lines(path):
    with (gzip.open(path,"rt") if path.suffix==".gz" else path.open()) as f:
        return [json.loads(s) for s in f if s.startswith("{")]


def state_key(path):
    with gzip.open(path,"rt") as f:
        d=json.load(f)
    parts=[d[k] for k in ("active_material_states","load","secant","displacement_m")]
    key=hashlib.sha256(json.dumps(parts,separators=(",",":")).encode()).hexdigest()
    assert key==d["state_hash"]
    return key


class Evidence(unittest.TestCase):
    def test_face_stop_only(self):
        for sign,n in (("plus",12),("minus",2)):
            old=read(P12/"c02"/(sign+"-summary.json"))
            r=read(OUT/"c02"/(sign+".json"))
            self.assertEqual(r["step"],old["consecutive"]+1)
            self.assertEqual(r["previous_radius_m"],old["stop"]["previous_radius_m"])
            self.assertEqual(r["threshold"],old["threshold"])
            self.assertEqual(r["radius_m"],1e-8)
            self.assertEqual(r["target_children"],n)
            kids=r["signed_children"]
            self.assertEqual(len({x["child"] for x in kids}),n)
            self.assertEqual({int(k) for k in r["child_directional_tangent_lower_integrals_kn_m"]},
                             {x["child"] for x in kids})
            for x in kids:
                self.assertEqual(x["leaves"],4**r["depth"])
                self.assertEqual(np.array(x["signed_nodal_center_kn"]).shape,(15,2))
                self.assertEqual(np.max(np.abs(x["signed_nodal_past_kn"])),0)
            self.assertEqual(r["unscaled_included"],r["lower_kn_m"]*r["radius_m"]>r["center_residual_interval_kn"][1]+r["paired_past_upper_kn"])
            inv=r["inverse_test"]
            with np.load(OUT/"c02"/(sign+"-residual.npz")) as vec:
                raw=np.linalg.norm(abs(vec["center_kn"])+vec["center_radius_kn"])
                applied=np.linalg.norm(abs(vec["after_inverse_m"])+vec["after_inverse_radius_m"])
                self.assertLess(raw,r["center_residual_interval_kn"][1])
                self.assertLess(applied,inv["center_after_inverse_upper_m"])
            self.assertEqual(inv["included"],inv["contraction_upper"]<1 and inv["ball_image_upper_m"]<r["radius_m"])
            self.assertEqual(inv["ball_image_upper_m"],inv["center_after_inverse_upper_m"]+r["radius_m"]*inv["contraction_upper"])
            self.assertEqual(r["included"],r["unscaled_included"] or inv["included"])
            self.assertEqual(r["consecutive"],old["consecutive"]+int(r["included"]))
            self.assertIsNone(r["reequilibrated_full_history_relative_upper"])

    def test_reload(self):
        folder=OUT/"c01/reload"
        rpt=read(folder/"comparison.json")
        rows=lines(folder/"trace.jsonl.gz")
        old=[r for r in lines(P12/"c01/reload/entrance.jsonl.gz") if r.get("stage")=="reload" and r["accepted"]]
        trials=[r for r in rows if r.get("stage")=="reload"]
        self.assertEqual(old,[r for r in trials if r["accepted"]][:len(old)])
        self.assertEqual(len(trials),rpt["budget"]["p13"]["attempts"])
        self.assertLessEqual(len(trials),1787)
        halves={(r["span"][0],r["span"][-1]):r for r in rows if "labels" in r}
        for s in rpt["switches"]:
            r=trials[s["attempt"]-1]
            self.assertTrue(r["accepted"])
            span=halves[r["from"],r["to"]]
            i=s["half"]-1
            self.assertEqual(span["span"][i:i+2],s["span"])
            self.assertEqual([v[30] for v in span["labels"][i:i+2]],s["labels"])
        for t,paths in rpt["common"].items():
            for label,file in (("ordinary",P11/"c01/reload/ordinary.jsonl.gz"),("p11",P11/"c01/reload/cross.jsonl.gz"),
                               ("p12",P12/"c01/reload/entrance.jsonl.gz"),("p13",folder/"trace.jsonl.gz")):
                seq=[r for r in lines(file) if r.get("stage")=="reload"]
                seq=seq[:paths[label]["attempts"]]
                self.assertEqual(seq[-1]["to"],float(t))
                self.assertTrue(seq[-1]["accepted"])
                self.assertEqual(sum(r["accepted"] for r in seq),paths[label]["accepted"])

    def test_knots_and_cost(self):
        folder=OUT/"c01/alternate"
        p12=read(P12/"c01/alternate/plan.json")
        plans=[read(folder/f"point{k}/plan.json") for k in (14,26)]
        a,b=(set(p["added"]) for p in plans)
        self.assertFalse(a & b)
        self.assertEqual(a|b,set(p12["added"]))
        for p in plans:
            self.assertEqual(set(p["knots"]),set(p12["p11_knots"])|set(p["added"]))
            for n in (256,512,1024):
                charges=read(folder/f"point{p['point']}"/str(n)/"charges.json")
                self.assertTrue(all(v["outer"]+v["extra"]==v["total"]<=1787 for v in charges.values()))

    def test_two_native_field_gates(self):
        for label in ("point14","point26"):
            folder=OUT/"c01/alternate"/label
            s=read(folder/"score.json")
            vals={}
            for n in (256,512,1024):
                sub=folder/str(n)
                old=read(P10/str(n)/"manifest.json")["fields"]
                meta=read(sub/"manifest.json")
                with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as p, zipfile.ZipFile(P10/str(n)/"prefix-xor.zip") as x, \
                     zipfile.ZipFile(P10/str(n)/"fields-after-open.zip") as t, zipfile.ZipFile(sub/"fields-xor-p10.zip") as d:
                    for name,m in old.items():
                        data=np.frombuffer(p10_bytes(n,m,name,p,x,t),"u1")
                        size=meta["count"]*int(np.prod(m["shape"]))*np.dtype(m["dtype"]).itemsize
                        data=data[:size]
                        if name+".bin" in d.namelist():
                            data=np.bitwise_xor(data,np.frombuffer(d.read(name+".bin"),"u1"))
                        self.assertEqual(hashlib.sha256(data).hexdigest(),meta["fields"][name]["sha256"])
                        vals.setdefault(name,[]).append(data.view(m["dtype"]).reshape((meta["count"],*m["shape"]))[:s["common_coordinates"]])
            passed=0
            for name,(a,b,c) in vals.items():
                ok=(np.array_equal(a,b) and np.array_equal(b,c)) if name in CATEGORIES else observed(a,b,c)["trend"] in ("decreasing_aligned","identical")
                self.assertEqual(ok,s["fields"][name]["full"]["passes"])
                passed+=ok
                if name=="plastic_slip":
                    dot=((a-b)*(b-c)).reshape(len(a),32,3).sum(-1)
                    with np.load(folder/"signed-32.npz") as f:
                        np.testing.assert_array_equal(dot,f[label])
            self.assertEqual((len(vals),passed),(22,s["passed_fields"]))
            self.assertEqual(s["passed"],passed==22 and s["common_coordinates"]==1792)

    def test_contributions(self):
        folder=OUT/"c01/alternate"
        rpt=read(folder/"comparison.json")
        for label in ("p11","p12","point14","point26"):
            src=P12/"c01/alternate" if label in ("p11","p12") else folder/label
            with np.load(src/"signed-32.npz") as f:
                dot=f[label][:rpt["common_coordinates"]]
            for j,stage in enumerate(STAGES):
                y=dot[j*256:(j+1)*256]
                if len(y):
                    s=rpt["stages"][stage][label]
                    self.assertEqual(y.sum(),s["sum_m2"])
                    np.testing.assert_array_equal(y[-1],s["endpoint"])
        matches=lines(folder/"common-recontact.jsonl.gz")
        self.assertEqual(len(matches),rpt["common_recontact"])
        for label in ("p11","p12","point14","point26"):
            if matches[-1]["t"]==1:
                np.testing.assert_array_equal(matches[-1][label],rpt["stages"]["recontact"][label]["endpoint"])

    def test_signed_half_steps(self):
        folder=OUT/"c01/alternate"
        report=read(folder/"intervals-7-9.json")
        self.assertEqual(report["points"],[7,9])
        grid=np.array(report["common"])
        sums={k:np.zeros((len(grid)-1,2)) for k in report["routes"]}
        for r in lines(folder/"signed-halves-7-9.jsonl.gz"):
            l,h=r["span"]
            j=grid.searchsorted(h,side="left")-1
            self.assertTrue(grid[j]<=l<h<=grid[j+1])
            v=(np.array(r["weight_m"])*r["delta_m"]).sum(-1)
            np.testing.assert_array_equal(v,r["term_m2"])
            sums[r["path"]][j]+=v
        endpoint=read(folder/"comparison.json")["stages"]["recontact"]
        for label,data in report["routes"].items():
            np.testing.assert_array_equal(sums[label],data["terms_m2"])
            np.testing.assert_allclose(data["initial_m2"]+sums[label].sum(0),data["endpoint_m2"],rtol=1e-8,atol=1e-28)
            if grid[-1]==1:
                np.testing.assert_allclose(np.array(endpoint[label]["endpoint"])[[7,9]],data["endpoint_m2"],rtol=1e-13,atol=1e-30)

    def test_nc_stopped_directions(self):
        folder=OUT/"nc"
        cmp=read(folder/"comparison.json")
        exact=read(folder/"direction-check.json")
        for label,entries in exact.items():
            self.assertTrue(entries)
            self.assertTrue(any(not e["passes"] for e in entries))
            for e in entries:
                self.assertEqual(e["passes"],abs(e["native_dx"]-e["exact_direction"])<=1e-14+5e-9*abs(e["native_dx"]))
            self.assertTrue(cmp["paths"][label]["movement"][-1]["unresolved"])

    def test_nc_record_before_solve(self):
        folder=OUT/"nc"
        arc=read(HERE/"design.json")["nc"]["common_arc_m"]
        rpt=read(folder/"comparison.json")
        rows=lines(folder/"global.jsonl")
        start=read(P12/"nc/comparison.json")["cumulative_including_replay"]
        for path in ("ordinary","predicted"):
            key=state_key(P12/"nc"/(path+"-end.json.gz"))
            with gzip.open(P12/"nc"/(path+"-end.json.gz"),"rt") as f:
                checkpoint=json.load(f)
            next_width=checkpoint["next_width"]
            lam=checkpoint["load"]
            seq=[r for r in rows if r["path"]==path]
            used=accepted=0
            width=None
            length=0.
            for r in seq:
                if "frozen_prediction" in r:
                    self.assertIsNone(width)
                    f=r["frozen_prediction"]
                    self.assertEqual(r["parent_hash"],key)
                    self.assertEqual((f["ordinary"],f["load"]),(next_width,lam))
                    w=min(f["ordinary"],.9*min([*f["partition_boundaries"].values(),*[v["left"]["width"] for v in f["sign_brackets"]]]))
                    unique=all(.9*b["left"]["width"]>=w for b in f["unresolved_brackets"])
                    self.assertEqual(f["unique_width"],unique)
                    self.assertEqual(f["predicted"],w if unique else None)
                    width=min(arc-length,f["ordinary"] if path=="ordinary" else w)
                elif "accepted" in r:
                    self.assertEqual((r["parent_hash"],r["width"]),(key,width))
                    used+=1
                    self.assertEqual(r["cumulative"],start+used)
                    if not r["accepted"]:
                        width/=2
                else:
                    accepted+=1
                    self.assertLess(r["load"],lam)
                    next_width,lam=r["next_width"],r["load"]
                    length+=width
                    self.assertEqual(length,r["arc"])
                    width=None
                    key=r["state_hash"]
            s=rpt["paths"][path]
            self.assertEqual((s["accepted"],s["attempts"],s["arc"]),(accepted,used,length))
            self.assertEqual(state_key(folder/(path+"-end.json.gz")),key)
            if s["stop"] is None:
                self.assertEqual(arc,length)
            start+=used
        self.assertEqual(start,rpt["cumulative"])
        self.assertLessEqual(start,2000)


if __name__=="__main__":
    unittest.main(verbosity=2)
