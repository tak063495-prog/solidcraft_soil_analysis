"""Exercise the C-03 independent constitutive return on the saved physical GP inputs."""

import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "interim/2026-12-18-c03-plastic-bulk-cause-diagnosis"
sys.path.insert(0, str(STUDY))
import c03diag_mc as mc


class EdgeClassification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        recorded = ROOT / "interim/2026-11-23-ac-scope-decisions/results/c03.json"
        cls.samples = json.loads(recorded.read_bytes())["points"]
        ground = json.loads((ROOT / "examples/ground-joint-water-mc.json").read_bytes())
        cls.material = dict(ground["materials"][0]["mohr_coulomb"], cohesion_kpa=30.0)

    def replay(self, row, resolved):
        state = row["previous_state"]
        trace = []
        strain = np.array(row["total_compressive_strain"]) / row["count"]
        for _ in range(row["count"]):
            state = mc.step(self.material, state, strain, trace=trace, resolve_edges=resolved)
        return state["last_return"], trace[-1]

    def test_same_input_on_both_sides_at_gp273(self):
        rows = [s for s in self.samples if s["gp"] == 273]
        self.assertEqual(len(rows), 6)
        for row in rows:
            with self.subTest(route=row["route"]):
                old, prior = self.replay(row, False)
                new, last = self.replay(row, True)
                self.assertEqual((old, new), ("shear_face", "extension_edge"))
                self.assertGreater(last["sweeps"], prior["sweeps"])
                for sweep in (last, prior):
                    self.assertEqual(sweep["threshold_kpa"], 1e-10 * sweep["scale_kpa"])
                if row["count"] == 1:
                    self.assertEqual(last["threshold_kpa"], prior["threshold_kpa"])
                self.assertLessEqual(
                    abs(last["principal_kpa"][0] - last["principal_kpa"][1]), last["threshold_kpa"]
                )

    def test_distinct_sides_of_yield_and_genuine_face_are_preserved(self):
        rows = [s for s in self.samples if s["gp"] in (152, 182) and s["route"].startswith("m1c1-")]
        self.assertEqual(len(rows), 4)
        for row in rows:
            with self.subTest(gp=row["gp"], route=row["route"]):
                self.assertEqual(self.replay(row, True)[0], row["label"])

    def test_active_minimum_and_apex_require_their_own_principal_gaps(self):
        # Pair order (major,minor): 01,02,10,12,20,21.
        compression = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        self.assertTrue(mc.active_edge_resolved(np.array([5.0, 1.0, 1.0]), compression, 1e-10))
        self.assertFalse(mc.active_edge_resolved(np.array([5.0, 1.01, 1.0]), compression, 1e-10))
        apex = np.array([1.0, 1.0, 0.0, 1.0, 0.0, 0.0])
        self.assertTrue(mc.active_edge_resolved(np.ones(3), apex, 1e-10))
        self.assertFalse(mc.active_edge_resolved(np.array([1.6e-10, 0.8e-10, 0.0]), apex, 1e-10))
        self.assertFalse(mc.active_edge_resolved(np.array([1.01, 1.0, 1.0]), apex, 1e-10))
