"""Verify event4 common-grid comparisons, route receipts and scoped evidence."""

import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parents[1] / "interim/2026-11-24-c01-event4-retry"
sys.path.insert(0, str(HERE))
from c01_retry_metrics import History, align, compare_pair, observed, sha256
from c01_retry_fields import CATEGORIES, STAGES


class Event4Retry(unittest.TestCase):
    def test_independent_friction_event_retries_and_finishes_the_cycle(self):
        from run_case import Methods, context

        study, mesh, ground, case, control = context(4)
        control.update(localize_friction_events=True, semismooth_contact=True, contact_history_bisection=True)
        result = study.reference._reference(
            mesh, ground, case, 4, 1, 7, False, True, maximum_m=5e-5,
            path_error=control, methods=Methods(None, None, True, None, None),
        )
        self.assertEqual((result["path"][-1]["stage"], result["path"][-1]["progress"]), ("return", 1))
        targets = [a for a in result["attempts"] if "friction_event_target" in a.get("joint_path_error", {})]
        self.assertTrue(targets)
        for a in targets:
            progress, gp = a["joint_path_error"]["friction_event_target"]
            self.assertFalse(a["accepted"])
            self.assertLess(a["progress_from"], progress)
            self.assertLess(progress, a["progress_to"])
            self.assertTrue(0 <= gp < 32)

    def test_signed_direction_and_zero_fields(self):
        values = np.arange(14, dtype=float).reshape(7, 2)
        linear = observed(values + 1, values + 0.5, values + 0.25)
        self.assertAlmostEqual(linear["observed_order"], 1)
        self.assertAlmostEqual(linear["cosine"], 1)
        self.assertEqual(observed(values + 1, values, values + 0.5)["trend"], "withheld")
        self.assertEqual(observed(values, values, values)["trend"], "identical")
        self.assertIsNone(observed(values + 1, values, values)["observed_order"])
        fields = values.copy()
        fields[3, 1] += 2.0
        pair = compare_pair(values, fields, 1)
        self.assertEqual(pair["maximum_location"], [STAGES[3], 1.0, 1])
        self.assertEqual(pair["absolute_l2"], 2.0)
        self.assertEqual(compare_pair(values, fields, 1, True)["mismatches"], 1)

    def test_alignment_of_partial_grid_and_corrupt_stream(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            data = folder / "u.bin"
            x = np.arange(1, 12, dtype="<f8").reshape(-1, 1)
            data.write_bytes(x.tobytes())
            cfg = dict(
                n=4,
                side="native",
                count=len(x),
                fields={"u": dict(shape=[1], dtype="<f8", sha256=sha256(data))},
            )
            (folder / "manifest.json").write_text(json.dumps(cfg), encoding="utf-8")
            stages = [dict(last_accepted=t) for t in (1, 1, 0.75, None, None, None, None)]
            history = History(folder, "native", 4, stages)
            np.testing.assert_array_equal(align([history], "u", 2)[0].ravel(), [2, 4, 6, 8, 10])
            with self.assertRaisesRegex(ValueError, "non-nested"):
                history.grid("u", 3)
            data.write_bytes(x[::-1].tobytes())
            with self.assertRaisesRegex(ValueError, "invalid field"):
                History(folder, "native", 4, stages)

    def test_fixed_scope_receipt(self):
        pin = json.loads((HERE / "ledger.json").read_bytes())
        subprocess.run(
            ["git", "merge-base", "--is-ancestor", pin["commit"], "HEAD"], cwd=HERE, check=True
        )
        target = pin["commit"] + ":" + pin["path"]
        blob = subprocess.check_output(["git", "rev-parse", target], cwd=HERE).decode().strip()
        raw = subprocess.check_output(["git", "cat-file", "blob", target], cwd=HERE)
        self.assertEqual(blob, pin["git_blob_sha1"])
        self.assertEqual(hashlib.sha256(raw).hexdigest(), pin["sha256"])
        self.assertEqual(raw, (HERE / "results/validation.json").read_bytes())

    def test_published_regression_and_scopes(self):
        folder = HERE / "results"
        state = json.loads((folder / "validation.json").read_bytes())
        self.assertEqual((state["default_limit"], state["stage_limit"]), (2000, 2065))
        for f, digest in state["sha256"].items():
            self.assertEqual(sha256(folder / f), digest, f)
        self.assertTrue(state["regression"])
        self.assertFalse(state["c01_complete"])
        routes = json.loads((folder / "routes.json").read_bytes())
        self.assertEqual(len(routes), 32)
        controls = [r for r in routes.values() if r["previous_replay"] is not None]
        self.assertEqual(len(controls), 10)
        for r in controls:
            self.assertTrue(all(r["previous_replay"].values()))
        for name, r in routes.items():
            self.assertEqual(r["completed_stages"], sum(s["last_accepted"] == 1 for s in r["stages"]))
            for s in r["stages"]:
                self.assertLessEqual(s["attempts"], state["stage_limit"])
                self.assertEqual(s["attempts"] - s["accepted"], sum(s["rejections"].values()))
            if r["mode"] in ("local", "local-ws", "local-active"):
                raw = gzip.decompress((folder / (name + ".jsonl.gz")).read_bytes())
                self.assertEqual(hashlib.sha256(raw).hexdigest(), r["accepted_history_sha256"])
                self.assertEqual(len(raw.splitlines()), r["accepted"])
        fields = json.loads((folder / "fields.json").read_bytes())["refinements"]
        for item in fields.values():
            archive = item["slip_differences"]
            raw = gzip.decompress((folder / archive["file"]).read_bytes())
            ab, bc = np.frombuffer(raw, dtype=archive["dtype"]).reshape(archive["shape"])
            norms = [np.linalg.norm(ab), np.linalg.norm(bc)]
            slip = item["fields"]["plastic_slip"]["common"]
            self.assertAlmostEqual(slip["observed_order"], float(np.log2(norms[0] / norms[1])))
            self.assertAlmostEqual(slip["cosine"], float(np.vdot(ab, bc) / (norms[0] * norms[1])))
            self.assertEqual(len(item["fields"]), 22)
            self.assertEqual(item["complete"], item["common_times"] == 7 * 256)
            for name, f in item["fields"].items():
                if name in CATEGORIES:
                    self.assertNotIn("observed_order", f["common"])
                elif f["common"]["trend"] == "decreasing_aligned":
                    self.assertGreater(min(f["common"]["cosine"], f["common"]["observed_order"]), 0)
