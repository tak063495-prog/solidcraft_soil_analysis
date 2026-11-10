"""CI must preserve every discovered test and module fixture."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import run_ci_python_tests as runner


class PythonCIPartition(unittest.TestCase):
    def inventory(self):
        return {'test_nc_mean_stress_evidence': ['mean.Case.a', 'mean.Case.b'],
                'test_nc_spatial_evidence': ['spatial.Case.a'],
                'test_a': ['a.Case.a', 'a.Case.b'], 'test_b': ['b.Case.a'],
                'test_c': ['c.Case.a']}

    def test_union_is_exact_with_no_test_repeated(self):
        modules = self.inventory()
        result = runner.partition(modules)
        expected = sorted(x for values in modules.values() for x in values)
        actual = [x for values in result.values() for x in values]
        self.assertEqual(sorted(actual), expected)
        self.assertEqual(len(actual), len(set(actual)))
        self.assertEqual(set(result), set(runner.GROUPS))

    def test_whole_modules_keep_class_and_module_fixtures_together(self):
        modules = self.inventory()
        result = runner.partition(modules)
        for names in modules.values():
            self.assertEqual(sum(set(names) <= set(values) for values in result.values()), 1)
        self.assertEqual(result['python-mean-evidence'], modules['test_nc_mean_stress_evidence'])
        self.assertEqual(result['python-spatial-evidence'], modules['test_nc_spatial_evidence'])

    def test_missing_or_duplicate_discovery_cannot_succeed(self):
        for change in ('missing', 'duplicate', 'empty'):
            modules = self.inventory()
            if change == 'missing':
                del modules['test_nc_mean_stress_evidence']
            elif change == 'duplicate':
                modules['test_b'] = modules['test_a']
            else:
                modules['test_b'] = []
            with self.assertRaises(ValueError):
                runner.partition(modules)

    def test_loader_failure_is_not_counted_as_a_discovered_test(self):
        suite = unittest.TestLoader().loadTestsFromName('module_that_does_not_exist_for_ci_inventory')
        with self.assertRaisesRegex(ValueError, 'discovery failed'):
            list(runner.flatten(suite))


if __name__ == '__main__':
    unittest.main()
