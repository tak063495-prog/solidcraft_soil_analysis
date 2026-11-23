"""Exercise both admissible construction inputs through the production CLI."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile


ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "interim/2026-11-23-impl-root-causes-a-nc"
CLI = ROOT / "target/release" / ("soil-fem-3d.exe" if os.name == "nt" else "soil-fem-3d")


class TestConstructionTransfer(unittest.TestCase):
    def run_case(self, name):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp) / name
            subprocess.run([sys.executable, STUDY / (name + ".py"), CLI, folder], check=True)
            return json.loads((folder / "results.json").read_bytes())

    def test_tip_crossing_carries_closed_form_nodal_force(self):
        result = self.run_case("a02")
        self.assertEqual((result["nodes"], result["gp"], result["accepted_birth"]), (81, 32, 4))
        self.assertEqual(result["normal_reference_kpa"], -52.5)
        for row in result["rows"]:
            self.assertLessEqual(row["maximum_component_error_kn"], 1e-10)
            self.assertGreater(row["minimum_cone_margin_kpa"], 0)

    def test_saturated_cut_transfers_head_and_wall_load(self):
        result = self.run_case("a04")
        self.assertEqual((result["states"], result["fields"]), (24, 49))
        self.assertLess(result["comparisons"]["exterior"]["maximum_absolute"], 1e-10)
        self.assertEqual(result["old_case_exit"], 1)

    def test_same_branch_restart_retains_all_accepted_states_and_trials(self):
        evidence = STUDY / "observations.zip"
        validation = json.loads((STUDY / "validation.json").read_bytes())
        self.assertEqual(hashlib.sha256(evidence.read_bytes()).hexdigest(), validation["archive_sha256"])
        with zipfile.ZipFile(evidence) as z:
            for size in (18, 144):
                for engine in ("native", "independent"):
                    prefix = f"a01/{size}/{engine}"
                    direct, first, resumed = [json.loads(gzip.decompress(z.read(prefix + suffix + ".gz")))
                                               for suffix in ("", "-4", "-resumed")]
                    self.assertEqual(first["points"] + resumed["points"], direct["points"])
                    self.assertEqual(first["trials"] + resumed["trials"], direct["trials"])
                    self.assertEqual(resumed["restart"], direct["restart"])
                    self.assertEqual(direct["stop"], "equilibrated candidate differs from accepted contact branch")
                    self.assertTrue(all(row["result_checkpoint"] == row["parent"]
                                        for row in direct["trials"] if not row["accepted"]))
