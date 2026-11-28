"""Original NC-01 inputs, solver-only variants and the eight-system decision."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'interim/2026-11-23-impl-nc01-residual-methods'
METHODS = ('repeated_refinement', 'extended_residual', 'ruiz', 'gmres_polish')


class NcResidualMethods(unittest.TestCase):
    def test_six_inputs_preserve_original_problem(self):
        with tempfile.TemporaryDirectory() as tmp:
            original, replay = (Path(tmp) / name for name in ('cloud', 'pc'))
            prepare = [sys.executable, str(REPORT / 'prepare.py')]
            subprocess.run([*prepare, original, '--cloud-original'], check=True, capture_output=True)
            subprocess.run([*prepare, replay, '--source', original], check=True, capture_output=True)
            self.assertEqual(sorted(p.name for p in original.iterdir()), sorted(p.name for p in replay.iterdir()))
            for p in original.iterdir():
                self.assertEqual(p.read_bytes(), (replay / p.name).read_bytes())
            case = json.loads((replay / 'stages.json').read_bytes())
            self.assertEqual(case['settings']['relative_linear_tolerance'], 1e-10)
            self.assertEqual(len(json.loads((replay / 'mesh.json').read_bytes())['tetrahedra10']), 96)
            for method in METHODS:
                variant = json.loads((replay / f'stages-{method}.json').read_bytes())
                self.assertEqual(variant['settings'].pop('dense_residual_method'), method)
                self.assertEqual(variant, case)
            reuse = subprocess.run([*prepare, replay, '--cloud-original'], capture_output=True)
            self.assertNotEqual(reuse.returncode, 0)

    def test_recorded_eligibility_uses_all_original_systems(self):
        evidence = json.loads((REPORT / 'measurements.json').read_bytes())
        self.assertEqual(evidence['linear_tolerance'], 1e-10)
        modes = {'baseline', 'refined', *METHODS}
        seen = set()
        passing = set(METHODS)
        for system in evidence['saved']:
            key = (system['file'], system['index'])
            self.assertNotIn(key, seen)
            seen.add(key)
            raw = (ROOT / key[0]).read_bytes().splitlines(keepends=True)[key[1]]
            self.assertEqual(hashlib.sha256(raw).hexdigest(), system['sha256'])
            self.assertEqual(set(system['methods']), modes)
            for mode, observation in system['methods'].items():
                self.assertEqual(observation['passed'], observation['residual'] <= 1e-10)
                if not observation['passed']:
                    passing.discard(mode)
        self.assertEqual(len(seen), 8)
        self.assertEqual(set(evidence['eligible_methods']), passing)
        self.assertEqual(set(evidence['paths']), modes)
        self.assertEqual(len({p['process']['binary_sha256'] for p in evidence['paths'].values()}), 1)
        for path in evidence['paths'].values():
            self.assertLessEqual(path['max_accepted_residual'], 1e-10)
            self.assertGreaterEqual(path['target'], path['last_accepted'])
