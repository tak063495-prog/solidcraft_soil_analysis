"""Replay the single 18 x 9 B-02 investigation with the production library."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parents[1] / "interim/2026-11-23-b02-a1-factor-study"
spec = importlib.util.spec_from_file_location("b02_a1_factor_analysis", HERE / "analyze.py")
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


class A1FactorStudy(unittest.TestCase):
    def test_live_single_axis_controls_and_unmet_root(self):
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory)
            study.build(raw)
            report = study.summarize(raw)

        self.assertFalse(report["a1_passed"])
        self.assertEqual(report["relative_width"], 0.001)
        modes = report["modes"]
        self.assertEqual(set(modes), set(study.MODES))
        for layer in modes["baseline"]["layers"]:
            self.assertFalse(layer["qx"]["passed"])
            self.assertFalse(layer["mxx"]["passed"])
            # Compare the reported precision of the two independent endpoint evaluations.
            self.assertEqual([round(x, 6) for x in layer["direct_edge_endpoint"]],
                             [round(x, 6) for x in layer["solved"]])
            tied = layer["exact_nodes_tying_QX"]
            self.assertEqual(round(tied["centroid"] + 500, 6),
                             round(tied["centroid_minus_continuum_closed_form"], 6))
            self.assertEqual([round(x, 6) for x in layer["exact_both"]], [-500., 5000.])
            self.assertLess(layer["swap_readout"][0], layer["solved"][0])
        for mode in ("gauss_edge", "full_clamp", "single_section", "si_units",
                     "duffy_shear", "membrane_direct"):
            for value, baseline in zip(modes[mode]["A1_per_layer"], modes["baseline"]["A1_per_layer"]):
                self.assertEqual(round(value, 3), round(baseline, 3), mode)
