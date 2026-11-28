"""Solver-independent checks of the C-01 method controls."""

import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

HERE = Path(__file__).resolve().parents[1] / "interim/2026-11-23-c01-method-improvements"
sys.path.insert(0, str(HERE))
from methods import STAGES, Seed, Sphere, active_tangent


class C01Methods(unittest.TestCase):
    def test_published_counts_cover_all_methods_with_original_budget(self):
        directory = HERE / "results"
        validation = json.loads((directory / "validation.json").read_bytes())
        for name, key in [("matrix.csv", "counts_sha256"), ("routes.json", "routes_sha256")]:
            self.assertEqual(
                hashlib.sha256((directory / name).read_bytes()).hexdigest(), validation[key]
            )
        with (directory / "matrix.csv").open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        modes = {"base", "w", "s", "a", "ws", "wa", "sa", "wsa"}
        self.assertEqual(len(rows), 32)
        self.assertEqual(
            {(r["solver"], r["N"], r["mode"]) for r in rows},
            {(s, n, m) for s in ("native", "independent") for n in ("512", "1024") for m in modes},
        )
        for row in rows:
            counts = [int(row[stage]) for stage in STAGES]
            self.assertEqual(sum(counts), int(row["total_attempts"]))
            self.assertTrue(all(0 <= c <= 2000 for c in counts))
            self.assertEqual(row["reload_attempts"], row["reload"])
            if row["stop_stage"]:
                self.assertEqual(row[row["stop_stage"]], "2000")
                self.assertLess(float(row["stop_progress"]), 1.0)
            else:
                self.assertEqual(float(row["reload_accepted"]), 1.0)
        self.assertEqual(validation["stage_attempt_budget"], 2000)
        self.assertFalse(validation["c01_complete"])

    def test_source_interpolates_displacement_across_stage_origins(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "seed.jsonl"
            with path.open("w") as stream:
                for stage in range(7):
                    for time in (0.5, 1.0):
                        print(json.dumps([stage, time, [[stage + time, 0.0, 0.0]]]), file=stream)
            seed = Seed(path, 1)
            np.testing.assert_array_equal(seed.interpolate("geostatic", 0.25), [0.25, 0, 0])
            np.testing.assert_array_equal(seed.interpolate("reload", 0.25), [5.25, 0, 0])
            path.write_text(path.read_text().splitlines()[0] + "\n")
            with self.assertRaisesRegex(ValueError, "incomplete"):
                Seed(path, 1)

    def test_arc_correction_meets_sphere_after_equilibrium_retrial(self):
        arc = Sphere(1.0)
        u = np.zeros(3)
        retry = arc.retry("shear", 0.0, 1.0, 1.0, u, np.full(3, 0.75), False)
        self.assertEqual(retry, 0.8)
        self.assertIsNone(arc.retry("shear", 0.0, retry, 1.0, u, np.full(3, 0.6), False))
        self.assertEqual((arc.calls, arc.corrections), (2, 1))
        self.assertIsNone(arc.retry("reload", 0.0, 1.0, 1.0, u, np.ones(3), True))

    def test_active_branches_include_reclosing_fraction_derivative(self):
        system = SimpleNamespace(
            normal=np.array([[0.0, 0.0, 1.0]]),
            kn=np.array([8.0]),
            kt=np.array([4.0]),
            mu=np.array([0.5]),
        )
        initial = np.zeros((1, 3))
        for old, target, expected in [
            (initial, [0.0, 0.0, 0.25], np.zeros((3, 3))),
            (initial, [0.125, 0.0, -1.0], np.diag([4.0, 4.0, 8.0])),
            (initial, [1.0, 0.0, -0.25], [[0, 0, -4], [0, 1, 0], [0, 0, 8]]),
            (np.array([[0.0, 0.0, 0.25]]), [0.25, 0.0, -0.25], [[2, 0, -1], [0, 2, 0], [0, 0, 8]]),
        ]:
            jac = active_tangent(system, old, initial, np.array([target]), np.array([1]))
            np.testing.assert_array_equal(jac[0], expected)
