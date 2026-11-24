"""Analytical controls for independent reference drivers."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


DIR = (
    Path(__file__).resolve().parents[1]
    / "interim/2026-11-24-independent-local-references"
)


def load_module(name):
    spec = importlib.util.spec_from_file_location(name, DIR / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestLocalReferences(unittest.TestCase):
    def test_wall_requires_both_axes_including_boundary(self):
        wall = load_module("assess_wall")
        import json

        root = DIR.parent
        spec = json.loads(
            (root / "2026-11-24-t9d-benchmarks/criteria.json").read_bytes()
        )
        grids = json.loads(
            (root / "2026-11-24-approved-four/kirsch-wall.json").read_bytes()
        )
        self.assertTrue(wall.assess(grids, spec)["passed"])
        target = grids[-1]["stages"][-1]["extrapolation"][1]
        target["wall_kpa"] = 15450.0
        self.assertTrue(wall.assess(grids, spec)["passed"])
        target["wall_kpa"] = 15450.001
        self.assertFalse(wall.assess(grids, spec)["passed"])

    def test_rotated_mc_face_and_edges_satisfy_consistency(self):
        mc = load_module("c03_point")
        material = dict(
            young_kpa=100000.0,
            poisson=0.3,
            cohesion_kpa=30.0,
            friction_deg=30.0,
            dilation_deg=5.0,
            flow="non_associated",
            tension_cutoff_kpa=None,
        )
        law = mc.PointReturn(material)
        angle = 0.37
        rot = np.array(
            [
                [np.cos(angle), 0, np.sin(angle)],
                [0, 1, 0],
                [-np.sin(angle), 0, np.cos(angle)],
            ]
        )
        old = np.repeat((100 * np.eye(3))[None, :, :], 3, axis=0)
        de = np.array(
            [
                rot @ np.diag(v) @ rot.T
                for v in (
                    (-0.002, 0, 0.002),
                    (-0.002, -0.002, 0.004),
                    (-0.004, 0.002, 0.002),
                )
            ]
        )
        stress, dep, labels = law.step(old, de)
        np.testing.assert_array_equal(labels, [1, 2, 3])
        np.testing.assert_allclose(
            np.linalg.eigvalsh(stress) @ law.normal, law.k, rtol=1e-12
        )
        elastic = de - dep
        expected = (
            old
            + law.g2 * elastic
            + law.lame * np.trace(elastic, axis1=1, axis2=2)[:, None, None] * np.eye(3)
        )
        np.testing.assert_allclose(stress, expected, atol=1e-10)

    def test_elastic_triaxial_strain_and_pressure(self):
        mcc = load_module("mcc_history")
        import json

        spec = json.loads(
            (DIR.parent / "2026-11-24-t9d-benchmarks/criteria.json").read_bytes()
        )["mcc"]
        vals = mcc.reference(spec, 40, 2)
        p, q, v, ea, el, ev = vals[0]
        self.assertAlmostEqual(p - q / 3, spec["p0_kpa"])
        self.assertAlmostEqual(q / (3 * 250), ea - ev / 3)
        self.assertAlmostEqual(ev, ea + 2 * el)
        v0 = spec["v_lambda"] - spec["lambda"] * np.log(40) + spec["kappa"] * np.log(8)
        self.assertAlmostEqual(v, v0 - spec["kappa"] * np.log(p / 5))

    def test_recorded_decisions_follow_history_extrema_and_point_labels(self):
        import gzip
        import json

        mcc = load_module("mcc_history")
        measures = json.loads((DIR / "mcc.json").read_bytes())
        archived = json.loads(
            gzip.decompress((DIR / "mcc-history.json.gz").read_bytes())
        )
        self.assertEqual(len(measures["cases"]), 8)
        for case, raw in zip(measures["cases"], archived["histories"], strict=True):
            pc = str(int(case["pc_kpa"]))
            diffs = mcc.differences(
                np.asarray(raw["fem"]), np.asarray(archived["fine_reference"][pc])
            )
            self.assertEqual(diffs, case["maxima"])
            self.assertEqual(
                case["passed"], all(d["relative"] < 0.02 for d in diffs.values())
            )
        mc = load_module("c03_point")
        summary = json.loads((DIR / "c03.json").read_bytes())
        native = json.loads(gzip.decompress((DIR / "c03-inputs.json.gz").read_bytes()))
        refs = json.loads(gzip.decompress((DIR / "c03-history.json.gz").read_bytes()))
        for route, points in summary["routes"].items():
            frames = native[route]["frames"]
            cumulative = np.cumsum(
                [r["fem_work_increment_kpa"] for r in frames], axis=0
            )
            error = mc.discrepancies(
                cumulative, np.array(refs[route]["reference_running_work_kpa"])
            )
            for j, point in enumerate(points):
                np.testing.assert_allclose(error[j], point["relative_max"])
                matched = sum(
                    frame["fem_label"][j] == mc.LABELS[ref[j]]
                    for frame, ref in zip(
                        frames, refs[route]["reference_label"], strict=True
                    )
                )
                self.assertEqual(point["matched"], matched)
        self.assertEqual(
            load_module("scope_rows").rows(),
            json.loads((DIR / "scope-rows.json").read_bytes()),
        )
