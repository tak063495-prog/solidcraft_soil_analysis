"""Independent checks of the p17 evidence and the unmodified acceptance gates."""
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

HERE=Path(__file__).resolve().parent
OUT=Path(os.environ.get("P17_RESULTS",HERE))
P16=HERE.with_name("2026-11-28-remaining-p16")


def read(path):
    return json.loads(path.read_text())


class Evidence(unittest.TestCase):
    def test_certified_prefix(self):
        for sign in ("plus","minus"):
            a=read(P16/f"c02/{sign}.json")
            b=read(OUT/f"c02/{sign}.json")
            stop=a["steps"][-1]["step"]
            self.assertEqual([r for r in b["steps"] if r["step"]<stop],a["steps"][:-1])
            self.assertEqual((b["threshold"],b["total"]),(.001,4156))
            self.assertEqual(b["consecutive"],b["steps"][-1]["step"]-1)
            self.assertIsNone(b["full_history_relative_upper"])

    def test_only_largest_leaf_split(self):
        for path in (OUT/"c02").glob("*-leaves.npz"):
            with np.load(path) as d:
                rows=d["tile_remainder_m"]
                win=int(d["selected_index"])
                self.assertEqual(rows.shape[1],16)
                self.assertEqual(win,np.argmax(rows.sum(1)))
                self.assertEqual(d["child_tile_remainder_m"].size,64)
                self.assertEqual(len(np.unique(d["block_leaf"],axis=0)),len(rows))
                before=d["prior_direct_remainder_m"]
                np.testing.assert_allclose(rows.sum(),before,rtol=2e-13,atol=0.)
                update=before-rows[win].sum()+d["child_tile_remainder_m"].sum()
                np.testing.assert_allclose(update,d["updated_direct_remainder_m"],rtol=3e-14,atol=0.)
        self.assertTrue(all((OUT/f"c02/{tag}-leaves.npz").exists() for tag in ("plus-833","minus-832")))

    def test_ball_images(self):
        for sign in ("plus","minus"):
            steps=read(OUT/f"c02/{sign}.json")["steps"]
            for r in steps:
                if "p15_nodal_past" not in r: continue
                tag=f"{sign}-{r['step']}"
                with np.load(OUT/f"c02/{tag}-pair.npz") as p,np.load(OUT/f"c02/{tag}-direct.npz") as d:
                    c=d["center_after_inverse_m"]
                    base=c+p["mid_action_m"]
                    e=min(p["remainder_m"],np.linalg.norm(p["component_radius_m"]))
                    projection=np.sqrt(np.linalg.norm(base)**2+e**2+2*(np.linalg.norm(c)*p["center_direction_support_m"]+
                        np.linalg.norm(p["mid_action_m"])*e))
                    v=r["inverse_test"]
                    eta=min(projection,np.linalg.norm(base)+p["remainder_m"],
                        np.linalg.norm(abs(base)+p["component_radius_m"]))*(1+1e-9)+r["center_only_inverse_test"]["direct_center_remainder_m"]+v["force_guard_m"]
                    self.assertEqual(eta+1e-8*v["contraction_upper"],v["ball_image_upper_m"])
                    self.assertEqual(v["included"],v["contraction_upper"]<1 and v["ball_image_upper_m"]<1e-8)

    def test_reload(self):
        from reload_score17 import run
        root=OUT/"c01/reload"
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            (d/"trace.jsonl.gz").symlink_to((root/"trace.jsonl.gz").resolve())
            run(d)
            self.assertEqual(read(root/"comparison.json"),read(d/"comparison.json"))
        self.assertEqual(read(root/"prefix-identity.json")["common_canonical_sha256"],
                         read(P16/"c01/reload/prefix-identity.json")["common_canonical_sha256"])

    def test_return(self):
        path=OUT/"c01/alternate"
        summary=read(path/"summary.json")
        for route in ("outer","middle"):
            with np.load(path/f"{route}-return.npz") as d,np.load(P16/f"c01/alternate/{route}-1792.npz") as orig:
                t=d["progress"]
                ix=np.searchsorted(t,np.arange(1,257)/256)
                self.assertEqual(len(t),summary["common_endpoints"])
                for level in (256,512,1024):
                    spans=d[f"half_{level}"]
                    np.testing.assert_array_equal(spans[:,1],t)
                    self.assertTrue((spans[1:,0]>=t[:-1]).all())
                for k in ("u","jump","reaction"):
                    g=d[k+"_growth"]
                    sq=d[k+"_squared_norms"]
                    np.testing.assert_allclose(g.sum((1,2))+d[k+"_residual"],sq[1]-sq[0],
                        atol=3e-14*sq.max(),rtol=1e-9)
                    for p in (1,25,21):
                        np.testing.assert_allclose(g[ix,p,:].sum(),orig[k+"_point_growth"][-256:,p].sum(),rtol=1e-9)
        self.assertEqual(read(path/"gates.json"),read(P16/"c01/alternate/gates.json"))
        self.assertTrue(all(r["trace_identical"] for r in read(path/"replay-cost.json").values()))

    def test_global_trials(self):
        a=read(OUT/"nc/summary.json")
        self.assertEqual((a["start_count"],a["cumulative"],a["rtol"],a["atol"]),(1366,1368,5e-9,1e-14))
        for path,r in a["paths"].items():
            self.assertFalse(r["accepted"])
            self.assertEqual(r["parent_hash"],r["end_hash"])
            self.assertEqual(r["lambda_before"],r["lambda_after"])
            self.assertEqual(r["candidate_count"],30)
            self.assertEqual(r["additional_global_trials"],1)
            ds=(OUT/f"nc/{path}-directions.jsonl").read_text().splitlines()
            for row in map(json.loads,ds):
                self.assertLessEqual(abs(row["schur"]-row["refined"][-1]),1e-14+5e-9*abs(row["refined"][-1]))


if __name__=="__main__":
    unittest.main(verbosity=2)
