import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
import nc01b_inputs as nc


class PhysicalExcavation(unittest.TestCase):
    def test_physical_supports_and_seed_balance(self):
        for n in (12, 96):
            for profile in nc.PROFILES:
                with self.subTest(n=n, profile=profile):
                    mesh, ground, case = nc.inputs(profile, n)
                    xyz = np.array(mesh['nodes'])
                    self.assertEqual(len(mesh['tetrahedra10']), n)
                    for side, axis in [('xmin', 0), ('xmax', 0), ('ymin', 1), ('ymax', 1)]:
                        bound = (xyz[:, axis].min() if side.endswith('min') else xyz[:, axis].max())
                        self.assertEqual(set(mesh['node_sets'][side]), set(np.flatnonzero(xyz[:, axis] == bound)))
                    fixed = {(bc['node_set'], bc['direction'].index(1.)) for bc in case['boundary']['constraints']}
                    self.assertEqual(('xmax', 0) in fixed, profile != 'pressure')
                    self.assertEqual(('bottom', 0) in fixed, profile != 'k0')
                    geo = ground['initialization'].get('geostatic', {})
                    self.assertEqual([p['face_set'] for p in geo.get('base_loads', {}).get('surface_loads', [])],
                                     ['xmax'] if profile == 'pressure' else [])
                    result = nc.initial_balance(mesh, ground, case)
                    self.assertAlmostEqual(result['support_reaction_kn'][2], 1373.625, places=8)
                    self.assertLess(result['relative_free_residual'], ground['initialization']['relative_equilibrium_tolerance'])

    def test_pressure_reaction_is_required_before_construction(self):
        mesh, ground, case = nc.inputs('pressure', 12)
        ground['initialization']['geostatic']['base_loads']['surface_loads'].clear()
        with self.assertRaisesRegex(ValueError, 'not in equilibrium'):
            nc.initial_balance(mesh, ground, case)

    def test_profile_and_size_must_be_selected(self):
        with self.assertRaises(ValueError):
            nc.inputs('NC-01', 12)
        with self.assertRaises(ValueError):
            nc.inputs('roller', 768)
