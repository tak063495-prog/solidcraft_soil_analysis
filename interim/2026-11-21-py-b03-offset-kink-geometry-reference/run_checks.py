"""Run explicit B-03 checks and write machine-readable counts/IDs (no skip success)."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import unittest

import numpy
import scipy


def ids(suite):
    result = []
    for child in suite:
        if isinstance(child, unittest.TestSuite):
            result.extend(ids(child))
        else:
            result.append(child.id())
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir')
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    here = Path(__file__).resolve().parent
    suite = unittest.defaultTestLoader.discover(str(here), pattern='checks.py', top_level_dir=str(here))
    collected = ids(suite)
    start = time.monotonic()
    with (output/'new-checks.log').open('w', encoding='utf-8') as stream:
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    report = {
        'scope': 'new independent B-03 checks only', 'engineering_accepted': None,
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'skipped': len(result.skipped), 'expected_failures': len(result.expectedFailures),
        'unexpected_successes': len(result.unexpectedSuccesses), 'collected_ids': collected,
        'duplicate_ids': len(collected)-len(set(collected)),
        'elapsed_seconds': time.monotonic()-start, 'python': platform.python_version(),
        'numpy': numpy.__version__, 'scipy': scipy.__version__,
        'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted(here.glob('*.py'))},
    }
    report['checks_successful'] = (result.wasSuccessful() and not result.skipped and
                                   not result.expectedFailures and not result.unexpectedSuccesses and
                                   result.testsRun == len(collected) and len(set(collected)) == len(collected))
    (output/'new-checks.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k not in ('collected_ids', 'source_sha256')}, ensure_ascii=False))
    return 0 if report['checks_successful'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
