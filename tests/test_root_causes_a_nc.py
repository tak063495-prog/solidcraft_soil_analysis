"""Exercise both admissible construction inputs through the production CLI."""
import gzip
import hashlib
import importlib.util
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
    def run_case(self, name, *options):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp) / name
            subprocess.run([sys.executable, STUDY / (name + ".py"), CLI, folder, *options], check=True)
            return json.loads((folder / "results.json").read_bytes())

    def test_tip_crossing_carries_closed_form_nodal_force(self):
        result = self.run_case("a02")
        self.assertEqual((result["nodes"], result["gp"], result["accepted_birth"]), (81, 32, 4))
        self.assertEqual(result["normal_reference_kpa"], -52.5)
        for row in result["rows"]:
            self.assertLessEqual(row["maximum_component_error_kn"], 1e-10)
            self.assertGreater(row["minimum_cone_margin_kpa"], 0)

    def test_saturated_cut_transfers_head_and_wall_load(self):
        result = self.run_case("a04", "--reference", "decimal40")
        self.assertEqual(result["reference"], "decimal40")
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


SCOPE_STUDY = ROOT / "interim/2026-11-23-ac-scope-decisions"


class ScopeDecisions(unittest.TestCase):
    def test_nonuniform_unloaded_sliding_pads_transfer_all_components(self):
        spec = importlib.util.spec_from_file_location(
            "ac_scope_shear_transfer", SCOPE_STUDY / "shear_transfer.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        suffix = ".exe" if os.name == "nt" else ""
        probe = ROOT / "target/release/examples" / ("joint_a02_projection_probe" + suffix)
        subprocess.run(
            ["cargo", "build", "--release", "--locked", "--example", "joint_a02_projection_probe"],
            cwd=ROOT,
            check=True,
        )
        with tempfile.TemporaryDirectory() as tmp:
            result = module.solve(probe, Path(tmp) / "transfer")
        self.assertEqual([p["reference_kpa"] for p in result["pads"]], [[13, 0, -60], [0, -8, -40]])
        self.assertEqual((result["nodes"], result["gp"]), (81, 32))
        self.assertLessEqual(result["maximum_nodal_component_error_kn"], 1e-10)
        self.assertGreater(result["minimum_gp_cone_margin_kpa"], 0)

    def test_corrected_arrivals_depart_from_old_branch_on_the_same_points(self):
        data = json.loads((SCOPE_STUDY / "results/a01.json").read_bytes())
        for size, first_gp in (("18", 4), ("144", 11)):
            pairs = list(zip(data[size]["native"], data[size]["independent"], strict=True))
            self.assertEqual(len(pairs), 10)
            for native, independent in pairs:
                self.assertEqual(native["ds"], independent["ds"])
                self.assertEqual(native["accepted"], independent["accepted"])
                self.assertEqual(native["corrected"]["events"], independent["corrected"]["events"])
                self.assertEqual(native["corrected"]["stored_label_mismatches"], [])
            self.assertEqual(pairs[0][0]["parent_progress"], 0.02)
            self.assertTrue(pairs[0][0]["accepted"])
            self.assertEqual(pairs[1][0]["corrected"]["changed"][0]["gp"], first_gp)
            self.assertEqual(pairs[-1][0]["corrected"]["changed"][0]["candidate"], "slip")

    def test_scope_ledger_binds_results_and_keeps_accuracy_decisions_pending(self):
        scoped = json.loads((SCOPE_STUDY / "scope-rows.json").read_bytes())
        rows = scoped["scope_rows"]
        self.assertEqual(rows["a01_dry_to_002"]["status"], "達成")
        self.assertEqual(rows["a02_original"]["status"], "検証対象外")
        for key in ("a01_open_0125", "a01_wet_ownership", "a01_instability",
                    "a04_unsaturated", "a04_richards", "c02_accuracy", "c03_accuracy"):
            self.assertEqual(rows[key]["status"], "未達")
        self.assertEqual(set(scoped["evidence_sha256"]), {"a01", "a02", "c03"})
        for case, digest in scoped["evidence_sha256"].items():
            evidence = SCOPE_STUDY / "results" / (case + ".json")
            self.assertEqual(hashlib.sha256(evidence.read_bytes()).hexdigest(), digest)
        c03 = json.loads((SCOPE_STUDY / "results/c03.json").read_bytes())
        self.assertEqual(len(c03["full_fields"]), 7)
        for pair in c03["full_fields"].values():
            self.assertEqual((pair["accepted"], len(pair["fields"])), (224, 23))
