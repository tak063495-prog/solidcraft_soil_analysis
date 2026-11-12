"""PC-3転記値の限定受入、拒否条件、追記専用の診断出力。pwsh不要。"""
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from nc_pc_gate import check_rate_arc, inspect_arc, source_manifest


def fixture():
    return (dict(equilibrium_samples=128, stop_reason=None,
                 maximum_bordered_true_linear_residual=5.4291880423103136e-11,
                 passed=False), [dict(step=i) for i in range(128)])


class RateArcGate(unittest.TestCase):
    def test_pc3_values_accepted_without_production_acceptance(self):
        report, samples = fixture()
        before = copy.deepcopy((report, samples))
        self.assertEqual(check_rate_arc(report, samples), [])
        self.assertEqual((report, samples), before)

    def test_127_samples_rejected(self):
        report, samples = fixture()
        self.assertTrue(check_rate_arc(report, samples[:127]))
        report['equilibrium_samples'] = 127
        self.assertTrue(check_rate_arc(report, samples[:127]))

    def test_residual_exceeded_rejected(self):
        report, samples = fixture()
        report['maximum_bordered_true_linear_residual'] = 1.1e-10
        self.assertTrue(check_rate_arc(report, samples))

    def test_residual_boundary_accepted(self):
        report, samples = fixture()
        report['maximum_bordered_true_linear_residual'] = 1e-10
        self.assertEqual(check_rate_arc(report, samples), [])

    def test_stop_reason_strings_rejected(self):
        report, samples = fixture()
        for value in ('null', '', 'bounded step count; inconclusive', 'failure'):
            with self.subTest(value=value):
                report['stop_reason'] = value
                self.assertTrue(check_rate_arc(report, samples))

    def test_empty_samples_rejected(self):
        report, _ = fixture()
        self.assertTrue(check_rate_arc(report, []))

    def test_each_missing_report_key_rejected(self):
        report, samples = fixture()
        for key in report:
            with self.subTest(key=key):
                value = report.copy()
                del value[key]
                self.assertIn(f'report.{key}: missing', check_rate_arc(value, samples))

    def test_wrong_shapes_rejected(self):
        report, samples = fixture()
        for value in (None, {}, {'samples': samples}, [samples], [[{}]] * 128, [None] * 128):
            with self.subTest(shape=type(value).__name__):
                self.assertTrue(check_rate_arc(report, value))
        self.assertTrue(check_rate_arc([], samples))

    def test_passed_must_be_boolean_false(self):
        report, samples = fixture()
        for value in (True, None, 0, 'false', [], {}):
            with self.subTest(value=value):
                report['passed'] = value
                self.assertTrue(check_rate_arc(report, samples))

    def test_count_must_be_integer(self):
        report, samples = fixture()
        for value in (None, '128', 128.0, True, 0, 129):
            with self.subTest(value=value):
                report['equilibrium_samples'] = value
                self.assertTrue(check_rate_arc(report, samples))

    def test_residual_type_and_finiteness(self):
        report, samples = fixture()
        for value in (None, '5.429e-11', False, float('nan'), float('inf'), -1., [], {}):
            with self.subTest(value=value):
                report['maximum_bordered_true_linear_residual'] = value
                self.assertTrue(check_rate_arc(report, samples))

    def test_missing_and_malformed_files_are_diagnostic_failures(self):
        with tempfile.TemporaryDirectory() as temp:
            arc = Path(temp)
            self.assertFalse(inspect_arc(arc)['accepted'])
            (arc / 'report.json').write_text('{invalid', encoding='utf-8')
            (arc / 'samples.json').write_text('[]', encoding='utf-8')
            self.assertFalse(inspect_arc(arc)['accepted'])
            for value in ('NaN', '1e999'):
                (arc / 'report.json').write_text('{"residual":' + value + '}', encoding='utf-8')
                self.assertFalse(inspect_arc(arc)['accepted'])

    def test_cli_rc_diagnostics_and_no_overwrite_or_source_change(self):
        with tempfile.TemporaryDirectory(prefix='PC gate 日本語 ') as temp:
            root = Path(temp)
            arc = root / 'rate-arc400'
            arc.mkdir()
            report, samples = fixture()
            for name, value in (('report', report), ('samples', samples)):
                (arc / f'{name}.json').write_text(json.dumps(value), encoding='utf-8-sig')
            (root / 'pc3-failure.json').write_text('retained failure', encoding='utf-8')
            before = source_manifest(root)
            output = root / 'gate.json'
            command = [sys.executable, str(ROOT / 'tools/nc_pc_gate.py'), '--arc', str(arc), '--output', str(output)]
            child = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(child.returncode, 0, child.stderr)
            result = json.loads(output.read_text(encoding='utf-8'))
            self.assertTrue(result['accepted'])
            for key in ('passed', 'nc01_complete', 'p6_complete'):
                self.assertIs(result[key], False)
            saved = output.read_bytes()
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
            self.assertEqual(output.read_bytes(), saved)
            self.assertEqual(source_manifest(root), before)
            (arc / 'samples.json').write_text('[]', encoding='utf-8')
            rejected = root / 'rejected.json'
            command[-1] = str(rejected)
            self.assertEqual(subprocess.run(command, capture_output=True).returncode, 1)
            self.assertIn('samples:', json.loads(rejected.read_text(encoding='utf-8'))['failures'][0])
            self.assertEqual((root / 'pc3-failure.json').read_text(), 'retained failure')

    def test_reuse_manifest_detects_source_change_and_retains_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'rate-arc400').mkdir()
            source = root / 'rate-arc400/samples.json'
            source.write_text('[]', encoding='utf-8')
            failure = root / 'pc3-failure.json'
            failure.write_text('old failure', encoding='utf-8')
            saved = root / 'reuse.json'
            saved.write_text(json.dumps(dict(read_only_source=str(root), source_manifest=source_manifest(root))), encoding='utf-8')
            command = [sys.executable, str(ROOT / 'tools/nc_pc_gate.py'), '--verify-reuse', str(saved), '--output']
            self.assertEqual(subprocess.run(command + [str(root / 'same.json')], capture_output=True).returncode, 0)
            source.write_text('[{}]', encoding='utf-8')
            self.assertEqual(subprocess.run(command + [str(root / 'changed.json')], capture_output=True).returncode, 1)
            self.assertEqual(failure.read_text(), 'old failure')

    def test_powershell_gate_when_available(self):
        shell = os.environ.get('NC_GATE_POWERSHELL') or shutil.which('pwsh') or shutil.which('powershell')
        if not shell:
            self.skipTest('PowerShell unavailable; all Python gate checks still run')
        with tempfile.TemporaryDirectory(prefix='PC gate 日本語 ') as temp:
            command = [shell, '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File',
                       str(ROOT / 'interim/2026-12-08-pc3-rate-arc-gate-fix/test-rate-arc-gate.ps1'),
                       '-Python', sys.executable, '-Output', temp]
            child = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, encoding='utf-8')
            self.assertEqual(child.returncode, 0, child.stdout + child.stderr)
            result = json.loads((Path(temp) / 'shell-result.json').read_text(encoding='utf-8-sig'))
            self.assertEqual(result['failed'], 0)
            self.assertEqual(result['passed'], 10)


if __name__ == '__main__':
    unittest.main()
