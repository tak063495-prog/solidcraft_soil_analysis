"""Numerical schedule endpoints, monotonicity and exact half-waypoint."""
import importlib.util
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('nc_oblique_release_test', ROOT/'interim/2026-11-10-nc-oblique-and-joint-history/oblique.py')
oblique = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(oblique)


class ObliqueRelease(unittest.TestCase):
    def test_original_endpoints_and_waypoints(self):
        for waypoint in ((.8, .2), (.2, .8)):
            self.assertEqual(oblique.fractions(0., waypoint), (0., 0.))
            self.assertEqual(oblique.fractions(.5, waypoint), waypoint)
            self.assertEqual(oblique.fractions(1., waypoint), (1., 1.))

    def test_monotone_component_bounded_withdrawal(self):
        for waypoint in ((.8, .2), (.2, .8)):
            values = np.array([oblique.fractions(s, waypoint) for s in np.linspace(0, 1, 257)])
            self.assertTrue(np.all((0 <= values) & (values <= 1)))
            self.assertTrue(np.all(np.diff(values, axis=0) >= 0))

    def test_outside_original_component_box_is_rejected(self):
        for s, waypoint in ((-.01, (.8, .2)), (1.01, (.8, .2)), (.5, (1.1, .2)), (.5, (.8, -.1))):
            with self.assertRaises(ValueError):
                oblique.fractions(s, waypoint)
        # The historical Newton driver resolves output paths internally. Its
        # caller must use the same absolute root when closing recorded streams.
        with tempfile.TemporaryDirectory(dir=ROOT) as temp:
            args = SimpleNamespace(output=Path(os.path.relpath(Path(temp)/'result')))
            fake = SimpleNamespace(write=lambda path, value: None)
            def execute(controls):
                controls.output.mkdir(parents=True)
                fake.append(controls.output.resolve()/'samples.jsonl', dict(recording_smoke_test=True))
            fake.run = execute
            with patch.dict(sys.modules, {'bounded_study': fake}):
                oblique.run(args)
            self.assertTrue(args.output.is_absolute())
            self.assertTrue((args.output/'streaming-digests.json').is_file())
