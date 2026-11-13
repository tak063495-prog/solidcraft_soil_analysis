"""Run every discovered unittest in exactly one Windows job, by whole module.

The two NC evidence modules and three C-01 modules each get a job. Other modules
are sorted and round-robin partitioned, including tools tests. Assertions and
tolerances are unchanged; every job records the complete inventory and partitions.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import unittest

GROUPS = ('python-core-0', 'python-core-1', 'python-mean-evidence', 'python-spatial-evidence',
          'python-c01-switch', 'python-c01-path', 'python-c01-cause')
HEAVY = {'test_nc_mean_stress_evidence': 'python-mean-evidence',
         'test_nc_spatial_evidence': 'python-spatial-evidence',
         'test_c01_global_step_control': 'python-c01-switch',
         'test_c01_path_error_control': 'python-c01-path',
         'test_c01_cause_diagnosis': 'python-c01-cause'}


def partition(modules):
    if not modules or not set(HEAVY) <= set(modules):
        raise ValueError('complete discovered inventory must include all dedicated evidence modules')
    ids = [test for values in modules.values() for test in values]
    if not ids or len(ids) != len(set(ids)) or any(not values for values in modules.values()):
        raise ValueError('nonempty unique discovered tests required')
    shards = {group: [] for group in GROUPS}
    regular = sorted(set(modules)-set(HEAVY))
    for index, module in enumerate(regular):
        shards[f'python-core-{index % 2}'].extend(sorted(modules[module]))
    for module, group in HEAVY.items():
        shards[group] = sorted(modules[module])
    combined = [test for values in shards.values() for test in values]
    if sorted(combined) != sorted(ids) or len(combined) != len(set(combined)) or any(not x for x in shards.values()):
        raise ValueError('partition lost, repeated or emptied tests')
    return shards


def flatten(suite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from flatten(item)
        else:
            if isinstance(item, unittest.loader._FailedTest):
                raise ValueError('unittest discovery failed: '+item.id())
            yield item


def discover(root):
    tests = []
    for folder in ('tests', 'tools'):
        suite = unittest.TestLoader().discover(str(root/folder), pattern='test_*.py')
        tests.extend(flatten(suite))
    by_id = {}
    modules = {}
    for test in tests:
        name = test.id()
        if name in by_id:
            raise ValueError('duplicate discovered test: '+name)
        by_id[name] = test
        modules.setdefault(test.__class__.__module__, []).append(name)
    return by_id, partition(modules)


class TimedResult(unittest.TextTestResult):
    """Record wall time per test, including its own setup/teardown, for CI budgets."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.test_elapsed_seconds = {}

    def startTest(self, test):
        self.started = time.monotonic()
        super().startTest(test)

    def stopTest(self, test):
        self.test_elapsed_seconds[test.id()] = time.monotonic()-self.started
        super().stopTest(test)


class Tee:
    def __init__(self, log):
        self.log = log

    def write(self, value):
        sys.stdout.write(value)
        self.log.write(value)

    def flush(self):
        sys.stdout.flush()
        self.log.flush()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--group', choices=GROUPS, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--inventory-only', action='store_true')
    args = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    report = dict(group=args.group, executed=not args.inventory_only, exit_code=1)
    status = 1
    try:
        tests, shards = discover(root)
        selected = shards[args.group]
        report.update(discovered_test_ids=sorted(tests), all_group_test_ids=shards,
                      discovered_total=len(tests), assigned_total=len(selected),
                      all_discovered_tests_partitioned_exactly_once=True)
        if args.inventory_only:
            status = 0
        else:
            with (args.output/'python-tests.log').open('w', encoding='utf-8') as log:
                result = unittest.TextTestRunner(stream=Tee(log), verbosity=2, resultclass=TimedResult).run(unittest.TestSuite(tests[x] for x in selected))
            report.update(tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
                          skipped=len(result.skipped), expected_failures=len(result.expectedFailures),
                          unexpected_successes=len(result.unexpectedSuccesses),
                          test_elapsed_seconds=result.test_elapsed_seconds)
            status = int(not result.wasSuccessful() or result.testsRun != len(selected) or bool(result.skipped)
                         or bool(result.expectedFailures) or bool(result.unexpectedSuccesses))
    finally:
        report.update(exit_code=status, elapsed_seconds=time.monotonic()-started)
        (args.output/'python-test-inventory.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return status


if __name__ == '__main__':
    sys.exit(main())
