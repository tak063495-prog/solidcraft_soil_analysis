"""Compare the independent opt-in tangent to the actual contact projection."""

import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from verify_joint_coupling import contact

METHODS = ROOT / "interim/2026-11-23-c01-method-improvements/methods.py"
SPEC = importlib.util.spec_from_file_location("c01_recontact_methods", METHODS)
methods = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(methods)


def system(mu):
    return SimpleNamespace(
        normal=np.array([[0.0, 0.0, 1.0]]),
        kn=np.array([8.0]),
        kt=np.array([4.0]),
        mu=np.array([mu]),
    )


def traction(system, old, target, count):
    return contact(system, old, np.zeros_like(old), target, count)[0][0, 0]


def tangent(system, old, target, count):
    return methods.active_tangent(system, old, np.zeros_like(old), target, np.array([count]))[0]


def centered(system, old, target, count, h):
    return np.column_stack([
        (traction(system, old, target + offset, count)
         - traction(system, old, target - offset, count)) / (2 * h)
        for offset in h * np.eye(3)
    ])


class C01Recontact(unittest.TestCase):
    def check_exact_recontact(self, shear, closed_x):
        point = system(0.5)
        old = np.array([[0.0, 0.0, 0.00025]])
        target = np.array([[shear, 0.0, 0.0]])
        h = 1e-11
        dz = h * point.normal
        expected = np.zeros((3, 3))
        expected[:, 2] = [closed_x, 0.0, 8.0]
        for count in (1, 4):
            with self.subTest(count=count):
                at = traction(point, old, target, count)
                closed = (at - traction(point, old, target - dz, count)) / h
                opened = (traction(point, old, target + dz, count) - at) / h
                np.testing.assert_array_equal(at, np.zeros(3))
                np.testing.assert_array_equal(opened, np.zeros(3))
                np.testing.assert_allclose(closed, expected[:, 2], rtol=0, atol=1e-6)
                actual = tangent(point, old, target, count)
                np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
                np.testing.assert_allclose(actual[:, 2], closed, rtol=0, atol=1e-6)

    def test_exact_recontact_sliding(self):
        self.check_exact_recontact(0.0005, -4.0)

    def test_exact_recontact_sticking(self):
        self.check_exact_recontact(0.000125, -2.0)

    def test_frictionless_contact(self):
        point = system(0.0)
        initial = np.zeros((1, 3))
        old = np.array([[0.0, 0.0, 0.25]])
        expected = np.diag([0.0, 0.0, 8.0])
        for count in (1, 4):
            for shear in (0.0, 1e-6):
                with self.subTest(count=count, shear=shear):
                    target = np.array([[shear, 0.0, -1.0]])
                    np.testing.assert_allclose(
                        centered(point, initial, target, count, 1e-6), expected, rtol=0, atol=1e-8
                    )
                    np.testing.assert_array_equal(tangent(point, initial, target, count), expected)
                    target[0, 2] = 0.0
                    np.testing.assert_array_equal(tangent(point, old, target, count), expected)
                    target[0, 2] = 0.125
                    np.testing.assert_array_equal(tangent(point, old, target, count), np.zeros((3, 3)))
            np.testing.assert_array_equal(
                tangent(system(0.5), initial, np.array([[0.0, 0.0, -1.0]]), count),
                np.diag([4.0, 4.0, 8.0]),
            )

    def test_smooth_recontact_matches_centered_difference(self):
        point = system(0.5)
        old = np.array([[0.0, 0.0, 0.25]])
        for shear in (0.125, 0.75):
            target = np.array([[shear, 0.0, -0.3125]])
            for count in (1, 4):
                with self.subTest(count=count, shear=shear):
                    np.testing.assert_allclose(
                        tangent(point, old, target, count),
                        centered(point, old, target, count, 1e-7),
                        rtol=0,
                        atol=1e-7,
                    )


if __name__ == "__main__":
    unittest.main()
