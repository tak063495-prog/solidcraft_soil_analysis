"""Guard CI coverage as Cargo auto-discovers additional integration targets."""
import json
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from run_ci_rust_tests import partition, commands


class Coverage(unittest.TestCase):
    def test_real_metadata_is_covered_once_and_new_targets_are_not_lost(self):
        meta = json.loads(subprocess.check_output(
            ['cargo', 'metadata', '--locked', '--no-deps', '--format-version', '1'], cwd=ROOT))
        package, = [p for p in meta['packages'] if Path(p['manifest_path']).resolve() == ROOT/'Cargo.toml']
        targets = package['targets']
        for inventory in [targets, targets+[dict(name='newly_added_test', kind=['test'])]]:
            groups = partition(inventory)
            names = [name for group in groups.values() for name in group]
            self.assertCountEqual(names, [t['name'] for t in inventory if 'test' in t['kind']])
            self.assertEqual(len(names), len(set(names)))
            for group, assigned in groups.items():
                cmd, = commands(group, groups)
                self.assertEqual([cmd[i+1] for i, arg in enumerate(cmd) if arg == '--test'], assigned)
                self.assertIn('--release', cmd)
                self.assertNotIn('--ignored', cmd)
                self.assertNotIn('--skip', cmd)

    def test_core_retains_unit_binary_doc_and_example_checks(self):
        core = commands('rust-core', {})
        self.assertEqual(core, [
            ['cargo', 'test', '--locked', '--release', '--lib', '--bins'],
            ['cargo', 'test', '--locked', '--release', '--doc'],
            ['cargo', 'build', '--locked', '--release', '--examples']])

    def test_ambiguous_or_empty_metadata_fails_closed(self):
        for targets in [[], [dict(name='x', kind=['test'])]*2]:
            with self.assertRaises(ValueError): partition(targets)
        with self.assertRaises(ValueError): commands('rust-integration-0', {'rust-integration-0': []})


if __name__ == '__main__':
    unittest.main()
