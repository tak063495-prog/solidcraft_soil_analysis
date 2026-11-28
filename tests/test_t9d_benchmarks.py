"""Scope decisions, analytic signs, and native excavation handoff for T9(d)."""
import gzip
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/"interim/2026-11-24-t9d-benchmarks"


def module(name):
    spec = importlib.util.spec_from_file_location("t9d_"+name, HERE/(name+".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class T9Benchmarks(unittest.TestCase):
    def test_registered_rows_reproduce(self):
        stored = json.loads((HERE/"ledger.json").read_bytes())
        self.assertEqual(module("assess").assess(HERE), stored)
        rows = {r["id"]:r for r in stored["rows"]}
        self.assertEqual(len(rows),len(stored["rows"]))
        self.assertFalse(rows["NC-01:overall"]["passed"])
        for name in ("arch:full-curve","kirsch:field","mcc:strain-curve","mcc:history-accuracy"):
            self.assertFalse(rows["T9d:"+name]["passed"])

    def test_arch_rejections_remain_in_the_record(self):
        curves = json.loads((HERE/"arch-results.json").read_bytes())
        trials = json.loads(gzip.decompress((HERE/"arch-results-trials.json.gz").read_bytes()))
        for n, methods in trials["cases"].items():
            for method,history in methods.items():
                rows = [dict(zip(trials["columns"],v,strict=True)) for v in history]
                accepted = [r for r in rows if r["accepted"]]
                self.assertEqual(len(accepted),curves[n][method]["accepted"])
                self.assertEqual(len(rows)-len(accepted),curves[n][method]["rejected"])
                self.assertTrue(all(r["true_linear_residual"] is None or
                                    r["true_linear_residual"]<=1e-10 for r in accepted))
        self.assertFalse(curves["12"]["beam"]["figure7_peak_and_reversal_passed"])
        self.assertTrue(curves["60"]["beam"]["figure7_peak_and_reversal_passed"])
        self.assertIsNone(curves["96"]["beam"]["fold_load"])

    def test_kirsch_wall_axes(self):
        k = module("kirsch")
        for axis, hoop in ((0,75000.),(1,15000.)):
            radius = np.eye(3)[axis]
            stress, change = k.solution(radius)
            self.assertEqual(stress[axis,axis],0.)
            self.assertEqual(stress[0,1],0.)
            self.assertEqual(stress[1-axis,1-axis],-hoop)
            self.assertLess(change[axis],0.)

    def test_excavation_transfers_core_traction(self):
        k = module("kirsch")
        solver = ROOT/"target/release"/("soil-fem-3d.exe" if os.name=="nt" else "soil-fem-3d")
        with tempfile.TemporaryDirectory() as directory:
            result = k.study(solver,Path(directory),1)
        start, half, end = result["stages"]
        self.assertEqual([start["release"],half["release"],end["release"]],[0,.5,1])
        self.assertGreater(np.linalg.norm(result["held_interface_kn"]),0)
        np.testing.assert_array_equal(end["interface_kn"],[0.,0.,0.])
        limit = k.SPEC["linear_relative_max"]
        for stage in [result["initial"],start,half,end]:
            self.assertLessEqual(stage["free_relative"],limit)
            self.assertLessEqual(stage["linear_relative"],limit)
        self.assertLess(start["max_displacement_change_m"],half["max_displacement_change_m"])
        self.assertLess(half["max_displacement_change_m"],end["max_displacement_change_m"])
