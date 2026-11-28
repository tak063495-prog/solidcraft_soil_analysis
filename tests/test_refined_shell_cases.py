"""Four explicitly scoped shell cases and the unresolved root section."""
import json
import os
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "interim/2026-11-23-b01-b02-refined-cases"
sys.path.insert(0, str(STUDY))
import cylinders
import refined


class RefinedShellCases(unittest.TestCase):
    def test_edf_connectivity_corner_restraints_and_load(self):
        grid = cylinders.grid()
        case, c = cylinders.cylinder_input("b01-edf364")
        self.assertEqual((len(grid["vertices"]), len(case["nodes"]), len(case["faces"])),
                         (209, 781, 364))
        edges = {}
        for orig, face in zip(grid["triangles"], case["faces"], strict=True):
            a, b, d, ab, bd, da = face["nodes"]
            self.assertEqual(set(orig), {a, b, d})
            for i, j, mid in ((a, b, ab), (b, d, bd), (d, a, da)):
                self.assertEqual(edges.setdefault(tuple(sorted((i, j))), mid), mid)
        self.assertEqual(len(set(edges.values())), 781 - 209)
        for name in refined.CYLINDERS:
            case, c = cylinders.cylinder_input(name)
            step = case["steps"][1]
            self.assertEqual(case["nodes"][c], [0, 0, 300])
            self.assertEqual(case["section"]["young_kpa"], 3000)
            self.assertEqual([(i, f) for i, f in enumerate(step["nodal_load"]) if f],
                             [(6 * c + 2, -0.00025)])
            free = set(step["free_dofs"])
            for i, (_, y, _) in enumerate(case["nodes"]):
                if y == 300:
                    self.assertEqual({j for j in range(6) if 6*i+j not in free}, {0, 2, 4})

    def test_live_cases_and_unmet_root(self):
        release = ROOT / "target/release"
        suffix = ".exe" if os.name == "nt" else ""
        if not (release / ("soil-fem-3d" + suffix)).exists() or not (
            release / "examples" / ("t10_layered_cantilever" + suffix)
        ).exists():
            subprocess.run(["cargo", "build", "--locked", "--release", "--bin", "soil-fem-3d",
                            "--example", "t10_layered_cantilever"], cwd=ROOT, check=True)
        live = refined.measure(release)
        names = {"b01-regular16", "b01-edf364", "b02-a3-12x6", "b02-a3-18x9", "b02-a1"}
        self.assertEqual(set(live), names)
        saved = {name: json.loads((STUDY / "results" / (name + ".json")).read_bytes())
                 for name in names}
        for name, count, probe, reference, width in [
            ("b01-regular16", 512, "c_uz_m", -1.8248e-5, 0.03),
            ("b01-edf364", 364, "c_uz_m", -1.8248e-5, 0.03),
            ("b02-a3-12x6", 288, "a3_dz_m", -3.90625e-5, 0.005),
            ("b02-a3-18x9", 648, "a3_dz_m", -3.90625e-5, 0.005),
        ]:
            for row in (saved[name], live[name]):
                self.assertEqual(row["elements"], count)
                gate = row[probe]
                self.assertEqual((gate["reference"], gate["relative_width"]), (reference, width))
                self.assertTrue(gate["passed"] and abs(gate["measured"] - reference) <= abs(reference)*width)
        for row in (live["b02-a1"], saved["b02-a1"]):
            self.assertEqual([g["grid"] for g in row["grids"]], ["12x6", "18x9"])
            for result in row["grids"]:
                self.assertEqual({r["offset_m"] for r in result["layers"]}, {-0.2, 0.2})
                for layer in result["layers"]:
                    for key, reference in (("qx_n_m", -500), ("mxx_nm_m", 5000)):
                        gate = layer[key]
                        self.assertEqual((gate["reference"], gate["relative_width"]), (reference, 0.001))
                        self.assertFalse(gate["passed"])
                        self.assertGreater(abs(gate["measured"] - reference), abs(reference)*0.001)
                root, = [s for s in result["sections"] if s["x_m"] == 0 and s["offset_m"] == 0.2]
                self.assertEqual(root["trace_from"], "right")
                self.assertEqual(root["points"][0]["qx_n_m"], result["layers"][1]["qx_n_m"]["measured"])
                self.assertEqual(len(result["sections"]), 10)
