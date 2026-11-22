"""Regression checks for the independent diagnostics, isolated from legacy module names."""
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "interim/2026-11-22-impl-missing-features-a01-a02-a04"


class MissingFeatures(unittest.TestCase):
    def run_study_python(self, source):
        subprocess.run([sys.executable, "-c", source], cwd=STUDY, check=True, timeout=120)

    def test_omitted_zero_and_shifted_richards_datum(self):
        self.run_study_python("""
from a04 import fixture, reference
import numpy as np
mesh, ground, case, dt = fixture()
case["joint_aperture_storage_m"] = 0.001
assert "head_datum_m" not in case["flow"]
base = reference(mesh, ground, case, dt, 2)
assert "head_datum_m" not in case["flow"]
for datum in (0., 1.5):
    case["flow"]["head_datum_m"] = datum
    shifted = reference(mesh, ground, case, dt, 2)
    for a,b in zip(base,shifted,strict=True):
        for field in ("heads","u","flux"):
            np.testing.assert_array_equal(a[field],b[field])
        assert a["volume_change"] == b["volume_change"]
""")

    def test_installed_branch_carries_the_first_birth_increment(self):
        self.run_study_python("""
from a01 import fixtures, Model, run
for choice, expected in (("centered", 0), ("installed_stick", 8)):
    model = Model(*fixtures.stacked(steps=8))
    rows = []
    try:
        run(model,8,lambda row, fields: rows.append(row),birth_predictor=choice)
    except ArithmeticError as err:
        assert "independent line search" in str(err)
    assert sum(r["stage"] == "birth" for r in rows) == expected
""")
