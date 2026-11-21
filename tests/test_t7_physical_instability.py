"""Six T7 decisions, via the same equilibrium/path routines used for A-01."""

import copy
import hashlib
import zipfile
import json
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np

HERE = (
    Path(__file__).resolve().parents[1]
    / "interim/2026-11-22-impl-t7-physical-instability"
)
sys.path.insert(0, str(HERE))
from continuation import (
    Controls,
    Evaluation,
    bifurcation,
    identity,
    simple_point,
    spectrum,
    trace,
)
from references import euler
from scan_accepted import scan


def polynomial(power):
    def evaluate(parent, z):
        x, lam = z
        return Evaluation(
            np.array([x**power - lam]),
            np.array([[power * x ** (power - 1)]]),
            np.array([-1.0]),
            dict(x=x, lam=lam),
            1.0,
        )

    return evaluate


class PhysicalInstability(unittest.TestCase):
    def test_1_iteration_limit_keeps_physical_d(self):
        result = trace(
            polynomial(2),
            {},
            [-1, 1],
            [1, -1],
            Controls(0.15, 1, 1e-12, 1e-10, 0, event_bisections=0),
        )
        self.assertEqual(
            (result["termination_reason"], result["physical_classification"]),
            ("Newton iteration limit", "D"),
        )
        self.assertFalse(result["trials"][0]["accepted"])
        self.assertEqual(result["last_accepted"], identity({}))

    def test_2_fold_and_monotone_displacement_control(self):
        ev = polynomial(2)
        c = Controls(0.08, 30, 1e-12, 1e-10, 10)
        for plane in (None, [1.0, 0.0]):
            r = trace(ev, {}, [-0.6, 0.36], [1, -1], c, local_control=plane)
            pts = np.array([p["z"] for p in r["accepted"]])
            self.assertLess(np.min(np.diff(pts[:, 1])), 0)
            self.assertGreater(np.max(np.diff(pts[:, 1])), 0)
            np.testing.assert_allclose(
                pts[:, 0] ** 2, pts[:, 1], rtol=0, atol=c.force_tol
            )

    def test_3_symmetric_single_double_and_indefinite_pivot(self):
        for k, neg in [
            (np.eye(2), 0),
            (-np.eye(2), 2),
            (np.array([[0, 1], [1, 0]]), 1),
        ]:
            a = spectrum(k, [1.0, 0.0], [0.0, 0.0])
            self.assertEqual((a["negative_eigen"], a["negative_ldlt"]), (neg, neg))
        self.assertEqual(
            spectrum(np.zeros((2, 2)), [1.0, 2.0], [0.0, 0.0])["nullity"], 2
        )

    def test_3_left_and_right_singular_spaces(self):
        s = spectrum(np.array([[0.0, 2.0], [0.0, 0.0]]), [4.0, 3.0], [0.0, 1.25])
        self.assertNotIn("negative_eigen", s)
        self.assertEqual(s["null_residual_norm"], 1.25)
        self.assertEqual(s["null_load_norm"], 3.0)
        np.testing.assert_allclose(s["left_weak_projector"], np.diag([0, 1]))
        np.testing.assert_allclose(s["right_weak_projector"], np.diag([1, 0]))

    def test_3_original_48_dof_stopped_residual(self):
        d = HERE.parent / "2027-01-04-a01-high-precision-nullspace"
        original = json.loads((d / "measurements.json").read_bytes())
        with zipfile.ZipFile(d / "evidence.zip.part000") as sealed:
            data = sealed.read("raw/double-with-both-svd.json")
            self.assertEqual(
                hashlib.sha256(data).hexdigest(),
                original["raw_files_sha256"]["double-with-both-svd.json"],
            )
            mats = json.loads(data)["records"]
            events = sealed.read(
                "source/formal/stacked__native__Pnative__Rnative__g8/native/progress.jsonl"
            )
        tangent = next(
            e
            for line in reversed(events.splitlines())
            if (e := json.loads(line))["kind"] == "tangent"
        )
        residual = np.asarray(tangent["reduced_residual"])
        spec = spectrum(
            np.asarray(mats["analytic"]["matrix"]), np.zeros_like(residual), residual
        )
        self.assertEqual((len(residual), spec["nullity"]), (48, 2))
        ref = original["matrices"]["analytic"]["true_null_projection_norms"][
            "saved_residual"
        ]
        self.assertAlmostEqual(spec["null_residual_norm"], float(ref), places=9)
        self.assertEqual(original["physical_instability_proved"], False)

    def test_4_simple_fold_degenerate_and_rejected_points(self):
        ev = polynomial(2)
        a = [-0.5, 0.25]
        b = [0.5, 0.25]
        self.assertEqual(simple_point(ev, {}, a, b, 1e-10)["kind"], "simple fold")
        self.assertEqual(simple_point(ev, {}, [-0.5, 1.0], b, 1e-10)["kind"], "D")
        self.assertEqual(
            simple_point(polynomial(3), {}, [-0.5, -0.125], [0.5, 0.125], 1e-10)[
                "kind"
            ],
            "D",
        )

        def event(parent, z):
            r = ev(parent, z)
            r.event = (z[0] > 0,)
            return r

        self.assertEqual(simple_point(event, {}, a, b, 1e-10)["kind"], "D")

    def test_4_pitchfork_two_distinct_branches(self):
        def ev(parent, z):
            x, p = z
            return Evaluation(
                np.array([x * (x * x - p)]),
                np.array([[3 * x * x - p]]),
                np.array([-x]),
                {},
                1.0,
            )

        result = bifurcation(
            ev, {}, [0.0, 0.0], [0.1, 0.05, 0.025], lambda p: np.zeros(1), 1e-11
        )
        self.assertEqual(result["kind"], "bifurcation")
        for sign, branch in zip((-1, 1), result["branches"]):
            self.assertTrue(
                all(sign * x > 0 and abs(x * x - p) < 1e-11 for x, p in branch)
            )

        samples = [dict(z=[0, load], state={}, event=()) for load in (-0.25, 0.25)]
        located = scan(ev, samples, 1e-11, lambda state: state)["branch_candidates"]
        self.assertEqual(len(located), 1)
        self.assertEqual(located[0]["kind"], "bifurcation")
        np.testing.assert_allclose(located[0]["z"], [0, 0], atol=1e-13)

    def test_5_restart_preserves_history_and_records_rejections(self):
        ev = polynomial(2)
        parent = {"generation": [1, 1]}
        old = copy.deepcopy(parent)
        args = (ev, parent, [-0.5, 0.25], [1, -1])
        opts = dict(ds=0.08, force_tol=1e-12, linear_tol=1e-10, iterations=12)
        full = trace(*args, Controls(steps=5, **opts))
        cut = trace(*args, Controls(steps=3, **opts))
        saved = cut["accepted"][-1]
        resumed = trace(
            ev,
            saved["state"],
            saved["z"],
            cut["trials"][-1]["direction"],
            Controls(steps=2, **opts),
            restart=json.loads(
                json.dumps(cut["restart"], default=lambda a: a.tolist())
            ),
        )
        np.testing.assert_array_equal(
            full["accepted"][-1]["z"], resumed["accepted"][-1]["z"]
        )
        bad = copy.deepcopy(cut["restart"])
        bad["payload"]["k"] *= 2
        with self.assertRaisesRegex(ValueError, "digest"):
            trace(
                ev,
                saved["state"],
                saved["z"],
                cut["trials"][-1]["direction"],
                Controls(steps=2, **opts),
                restart=bad,
            )
        self.assertEqual(parent, old)
        first = full["critical"][1]
        self.assertTrue(
            all("rhs" in it and "matrix" in it for it in first["iterations"])
        )
        self.assertEqual(len(full["trials"]), 5)

    def test_6_euler_mesh_and_mode(self):
        series = [euler(n) for n in (4, 8, 16)]
        self.assertTrue(
            all(
                a["relative_load_error"] > b["relative_load_error"]
                for a, b in zip(series, series[1:])
            )
        )
        self.assertLess(series[-1]["relative_load_error"], 2e-7)
        self.assertLess(series[-1]["mode_max_error"], 1e-7)
        self.assertAlmostEqual(series[-1]["mode"][-1], 1.0)

    def test_6_original_mesh_and_eight_children(self):
        # Each historical study imports its own fixture modules. Use the same
        # isolated entry point as the numerical runs, even in full discovery.
        cmd = (
            "import json; from parallel import fixtures, refine; "
            "mesh, _, _ = fixtures.stacked(steps=8); "
            "print(json.dumps([mesh, refine(mesh)]))"
        )
        mesh, fine = json.loads(
            subprocess.check_output([sys.executable, "-c", cmd], cwd=HERE)
        )
        self.assertEqual((len(mesh["tetrahedra"]), len(fine["tetrahedra"])), (18, 144))
        xyz = np.asarray(fine["nodes"])
        for i, cell in enumerate(mesh["tetrahedra"]):
            pts = xyz[cell]
            volume = abs(np.linalg.det(pts[1:] - pts[0]))
            for c in fine["tetrahedra"][8 * i : 8 * i + 8]:
                x = xyz[c]
                self.assertAlmostEqual(np.linalg.det(x[1:] - x[0]), volume / 8)
        for face in ("LowerSeam", "UpperSeam"):
            self.assertFalse(
                set(fine["node_sets"][face + ":a"])
                & set(fine["node_sets"][face + ":b"])
            )


if __name__ == "__main__":
    unittest.main()
