"""Input preservation and evidence binding for the opt-in dense stage solver."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'interim/2026-11-22-impl-linear-solver'
FIXTURES = ROOT / 'tests/fixtures/linear-refinement'


class DenseRefinement(unittest.TestCase):
    def test_pc_inputs_are_byte_preserved_and_only_solver_option_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, comparison = (Path(tmp) / name for name in ('source', 'comparison'))
            prepare = [sys.executable, str(REPORT / 'prepare_inputs.py')]
            subprocess.run([*prepare, str(source), '--cloud-original'], check=True, capture_output=True)
            for name in ('mesh', 'ground', 'stages'):
                p = source / (name + '.json')
                p.write_bytes(p.read_bytes().replace(b'\n', b'\r\n'))
            subprocess.run([*prepare, str(comparison), '--source', str(source)], check=True, capture_output=True)
            for name in ('mesh', 'ground', 'stages'):
                p = name + '.json'
                self.assertEqual((comparison / p).read_bytes(), (source / p).read_bytes())
            baseline = json.loads((source / 'stages.json').read_bytes())
            refined = json.loads((comparison / 'stages-refined.json').read_bytes())
            self.assertIs(refined['settings'].pop('dense_lu_refinement'), True)
            self.assertEqual(refined, baseline)

    def test_rejected_systems_and_corrected_solutions_match_the_sealed_measurement(self):
        report = json.loads((REPORT / 'measurements.json').read_bytes())
        solutions = (REPORT / 'corrected-solutions.jsonl').read_bytes()
        self.assertEqual(hashlib.sha256(solutions).hexdigest(), report['corrected_solutions_sha256'])
        by_system = {r['system_sha256']: r for r in map(json.loads, solutions.splitlines())}
        self.assertEqual(len(by_system), 8)
        for name, records in [('nc12.json', report['nc12']), ('nc96.jsonl', report['nc96_systems'])]:
            rows = (FIXTURES / name).read_bytes().splitlines(keepends=True)
            self.assertEqual(len(rows), len(records))
            for raw, record in zip(rows, records):
                key = hashlib.sha256(raw).hexdigest()
                self.assertEqual(key, record['sha256'])
                self.assertEqual(json.loads(raw)['native_true_relative_residual'], record['native_relative_residual'])
                self.assertEqual(by_system[key]['corrected_native_relative_residual'],
                                 record['corrected_native_relative_residual'])
                self.assertGreater(record['native_relative_residual'], 1e-10)
                self.assertLessEqual(record['corrected_native_relative_residual'], 1e-10)
                self.assertLessEqual(record['corrected_decimal80_relative_residual'], 1e-10)

        terminal = (REPORT / 'refined-terminal.jsonl').read_bytes()
        end = report['refined_terminal']['system']
        self.assertEqual(hashlib.sha256(terminal).hexdigest(), end['sha256'])
        self.assertGreater(end['native_relative_residual'], report['linear_tolerance'])
        self.assertEqual(json.loads(terminal)['native_true_relative_residual'], end['native_relative_residual'])
        self.assertFalse(report['refined']['complete'])
        self.assertFalse(report['nc01_complete'])


if __name__ == '__main__':
    unittest.main()
