"""Regression checks for the mesh budget and accepted dissipation comparisons."""

import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "interim/2026-11-23-c01-mesh-budget"
spec = importlib.util.spec_from_file_location("c01_budget_measure", HERE / "measure.py")
measure = importlib.util.module_from_spec(spec)
spec.loader.exec_module(measure)


class MeshBudget(unittest.TestCase):
    def test_independent_budget_counts_corner_cells_once(self):
        # The historical reference loads verification modules by their original short names.
        code = """
from run_case import context
study, mesh, _, _, _ = context(512)
limit = study.reference.stage_attempt_budget
assert limit(mesh, {}) == limit(mesh, {"attempt_budget_policy":"fixed"}) == 2000
cfg = {"attempt_budget_policy":"mesh_proportional"}
assert limit(mesh,cfg) == 2065
assert limit({"tetrahedra":[[0,1,2,3]]}, cfg) == 173
"""
        subprocess.run([sys.executable, "-c", code], cwd=measure.PREVIOUS, check=True)

    def test_order_uses_common_accepted_times_across_all_stages(self):
        h = {
            n: {(s, i / n): [s + i / n + 1 / n] for s in range(7) for i in range(1, n + 1)}
            for n in (256, 512, 1024)
        }
        result = measure.convergence(h)
        self.assertEqual(result["observed_order"], 1.0)
        self.assertEqual(result["common_times"], 7 * 256)
        self.assertEqual(result["adjacent"]["512-1024"]["components"], 7 * 512)
        del h[1024][6, 1.0]
        with self.assertRaisesRegex(ValueError, "missing exact"):
            measure.convergence(h)

    def test_published_counts_and_accepted_fields(self):
        path = HERE / "results"
        state = json.loads((path / "validation.json").read_bytes())
        for key in ("stages", "routes", "dissipation"):
            file = path / (key + (".csv" if key == "stages" else ".json"))
            self.assertEqual(hashlib.sha256(file.read_bytes()).hexdigest(), state[key + "_sha256"])
        routes = json.loads((path / "routes.json").read_bytes())
        orders = json.loads((path / "dissipation.json").read_bytes())
        for side in ("native", "independent"):
            for mode in ("base", "ws"):
                base = routes[f"{side}-512-{mode}"]
                fine = routes[f"{side}-1024-{mode}"]
                self.assertTrue(base["same_stage_records_as_fixed"])
                self.assertEqual(orders[side + "-" + mode] is not None, fine["complete"])
        for route in routes.values():
            raw = gzip.decompress((path / route["archive"]).read_bytes())
            self.assertEqual(hashlib.sha256(raw).hexdigest(), route["accepted_history_sha256"])
            self.assertEqual(len(raw.splitlines()), sum(s["accepted"] for s in route["stages"]))
            for s in route["stages"]:
                self.assertEqual(s["accepted"] + sum(s["rejections"].values()), s["attempts"])
                self.assertLessEqual(s["attempts"], state["stage_limit"])
        self.assertEqual((state["default_limit"], state["stage_limit"]), (2000, 2065))
        self.assertFalse(state["c01_complete"])
