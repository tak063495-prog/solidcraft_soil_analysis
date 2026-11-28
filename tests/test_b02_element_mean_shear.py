"""Opt-in shear projection, unchanged default acceptance, and checkpoint mode identity."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
STUDY_DIR = ROOT / "interim/2026-11-24-b02-a1-shear-mean"
spec = importlib.util.spec_from_file_location("b02_mean_study", STUDY_DIR / "study.py")
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)
from cylinders import standard


class ElementMeanShear(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["cargo", "build", "--release", "--locked", "--bin", "soil-fem-3d",
                        "--example", "t10_layered_cantilever"], cwd=ROOT, check=True)
        release = ROOT / "target/release"
        suffix = ".exe" if os.name == "nt" else ""
        cls.cli = release / ("soil-fem-3d" + suffix)
        cls.plate = release / "examples" / ("t10_layered_cantilever" + suffix)

    def test_fixed_grid_a1_and_a3(self):
        for grid in ("12x6", "18x9"):
            with self.subTest(grid=grid):
                pair = study.compare(self.plate, grid)
                for method, accepted in (("centroid", False), ("element_mean", True)):
                    observation = pair[method]
                    self.assertTrue(observation["a3"]["passed"])
                    self.assertEqual(observation["a3"]["relative_width"], .005)
                    self.assertEqual(len(observation["layers"]), 2)
                    for layer in observation["layers"]:
                        for name, ref in (("qx", -500), ("mxx", 5000)):
                            gate = layer[name]
                            self.assertEqual((gate["reference"], gate["relative_width"]), (ref, .001))
                            self.assertEqual(gate["passed"], accepted)
                        # Exact nodal kinematics must recover constant shear throughout
                        # the cell with the opt-in projection, including away from A1.
                        if method == "element_mean":
                            for q in layer["exact_nodes_gp_qx"]:
                                self.assertEqual(round(q, 6), -500.)
                        else:
                            self.assertGreater(layer["exact_nodes_center_qx"], -500.)

    def test_cli_checkpoint_is_bound_to_shear_choice(self):
        cards = json.loads(standard.CARDS.read_bytes())["benchmarks"]
        card = next(c for c in cards if c["id"] == "pinched_cylinder")
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            input_file = folder / "case.json"
            input_file.write_text(json.dumps(standard.build(card, 2, False)), encoding="utf-8")

            def run(name, args, check=True):
                return subprocess.run([self.cli, "curved-shell", input_file, folder/name, *args],
                                      capture_output=True, text=True, check=check)

            run("full", ["--shear-mean"])
            run("checkpoint", ["--shear-mean", "--until", "1"])
            run("resumed", ["--shear-mean", "--resume", folder/"checkpoint"])
            self.assertEqual((folder/"full").read_bytes(), (folder/"resumed").read_bytes())
            result = run("incompatible", ["--resume", folder/"checkpoint"], check=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("checkpoint identity", result.stderr)
