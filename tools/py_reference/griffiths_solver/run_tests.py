"""Run every new diagnostic test and reject skips or discovery omissions."""
import argparse
import json
from pathlib import Path
import sys
import time
import unittest


def flatten(suite):
    for value in suite:
        if isinstance(value, unittest.TestSuite):
            yield from flatten(value)
        else:
            yield value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists():
        parser.error('report already exists')
    suite = unittest.defaultTestLoader.discover(str(Path(__file__).resolve().parent), pattern='test_*.py')
    tests = list(flatten(suite))
    ids = [test.id() for test in tests]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError('nonempty unique discovery required')
    started = time.monotonic()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(tests))
    passed = bool(result.wasSuccessful() and result.testsRun == len(ids) and not result.skipped
                  and not result.expectedFailures and not result.unexpectedSuccesses)
    report = dict(discovered_test_ids=ids, discovered_total=len(ids), tests_run=result.testsRun,
                  failures=len(result.failures), errors=len(result.errors), skipped=len(result.skipped),
                  expected_failures=len(result.expectedFailures), unexpected_successes=len(result.unexpectedSuccesses),
                  elapsed_seconds=time.monotonic()-started, passed=passed)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return int(not passed)


if __name__ == '__main__':
    sys.exit(main())
