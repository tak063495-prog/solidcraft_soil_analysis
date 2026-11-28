"""実生成JSONの型とPC-4再利用の拒否条件。大規模経路は求解しない。"""
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT / 'interim/2026-11-09-nc-pseudo-arclength'
EVIDENCE = ROOT / 'interim/2026-12-05-nc-independent-arc-residual'
sys.path[:0] = [str(ROOT / 'tools'), str(PRIOR)]
import independent_arc as arc
import nc_pc_gate as gate

PC3_HEAD = '8342762bd714e33523ce5fef772572699fc21881'


class PcGateSchema(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        fixture = EVIDENCE / 'reference400-fixture.zip'
        expected = json.loads((EVIDENCE / 'reference400-fixture-sha256.json').read_text())
        if hashlib.sha256(fixture.read_bytes()).hexdigest() != expected[fixture.name]:
            raise AssertionError('reference fixture hash differs')
        with zipfile.ZipFile(fixture) as archive:
            archive.extractall(cls.root / 'reference400')
        # 生成器を本当に1点実行。128点の物理解・PC出力の再計算ではない。
        cls.report = arc.run(cls.root / 'reference400', cls.root / 'one-sample', 2e-5, 1)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def reuse_fixture(self, root):
        shutil.copytree(self.root / 'reference400', root / 'reference400')
        target = root / 'rate-arc400'
        shutil.copytree(self.root / 'one-sample', target)
        # 再利用検査専用の合成データ。同じ1点を複製したもので物理的受入に使わない。
        row = gate.read_json(target / 'samples.json')[0]
        rows = [dict(copy.deepcopy(row), step=i) for i in range(128)]
        (target / 'samples.json').write_text(json.dumps(rows), encoding='utf-8')
        report = dict(self.report, equilibrium_samples=128)
        (target / 'report.json').write_text(json.dumps(report), encoding='utf-8')
        controls = gate.read_json(target / 'controls.json')
        controls['count'] = 128
        (target / 'controls.json').write_text(json.dumps(controls), encoding='utf-8')
        (root / 'head.txt').write_text(PC3_HEAD, encoding='utf-8-sig')
        (root / 'pc3-failure.json').write_text('original failure retained', encoding='utf-8')
        for name, arguments in (
            ('reference400', ['tools/low_preload_nc_reference.py', '--steps', '400', '--method', 'rk4']),
            ('rate-arc400', ['interim/2026-11-09-nc-pseudo-arclength/independent_arc.py', '--count', '128', '--step-size', '0.00002']),
        ):
            (root / f'{name}-command.json').write_text(json.dumps(dict(head=PC3_HEAD, return_code=0, arguments=arguments)), encoding='utf-8-sig')
        return target

    def test_real_producer_writes_expected_json_names_and_types(self):
        report = gate.read_json(self.root / 'one-sample/report.json')
        rows = gate.read_json(self.root / 'one-sample/samples.json')
        self.assertIs(type(report['equilibrium_samples']), int)
        self.assertEqual(report['equilibrium_samples'], 1)
        self.assertIsNone(report['stop_reason'])
        self.assertIs(type(report['maximum_bordered_true_linear_residual']), float)
        self.assertLessEqual(report['maximum_bordered_true_linear_residual'], 1e-10)
        self.assertIs(report['passed'], False)
        self.assertIs(type(rows), list)
        self.assertEqual(len(rows), 1)
        self.assertIs(type(rows[0]), dict)
        self.assertTrue(gate.check_rate_arc(report, rows))  # 1点は128点ゲートを満たさない。

    def test_reuse_validation_keeps_source_bytes_and_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.reuse_fixture(root)
            before = gate.source_manifest(root)
            result = gate.inspect_reuse(root, ROOT)
            self.assertTrue(result['accepted'], result['failures'])
            self.assertIs(result['passed'], False)
            self.assertEqual(gate.source_manifest(root), before)
            self.assertEqual(result['source_manifest'], before)
            self.assertIn('pc3-failure.json', before)

    def test_reuse_refuses_bad_command_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.reuse_fixture(root)
            path = root / 'rate-arc400-command.json'
            original = gate.read_json(path)
            for changed in (dict(original, return_code=1), dict(original, return_code=False),
                            dict(original, head='wrong'), dict(original, arguments=[])):
                with self.subTest(changed=changed):
                    path.write_text(json.dumps(changed), encoding='utf-8')
                    self.assertFalse(gate.inspect_reuse(root, ROOT)['accepted'])

    def test_reuse_refuses_missing_reference_or_bad_seed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            target = self.reuse_fixture(root)
            rows = gate.read_json(target / 'samples.json')
            rows[0]['tangent'] = []
            (target / 'samples.json').write_text(json.dumps(rows), encoding='utf-8')
            self.assertFalse(gate.inspect_reuse(root, ROOT)['accepted'])
            (root / 'reference400/unload.json').unlink()
            self.assertFalse(gate.inspect_reuse(root, ROOT)['accepted'])

    def test_reuse_refuses_changed_source_code(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.reuse_fixture(root)
            check_output = subprocess.check_output
            def changed(command, **kwargs):
                data = check_output(command, **kwargs)
                if command[:2] == ['git', 'show'] and command[-1].endswith(':tools/mcc_rate_rk4.py'):
                    return data + b'\n# different source\n'
                return data
            with patch.object(gate.subprocess, 'check_output', changed):
                result = gate.inspect_reuse(root, ROOT)
            self.assertFalse(result['accepted'])
            self.assertTrue(any('source/input differs' in f for f in result['failures']))

    def test_real_rerun_saves_specific_failure_and_preserves_pc3(self):
        shell = os.environ.get('NC_GATE_POWERSHELL') or shutil.which('pwsh') or shutil.which('powershell')
        if not shell:
            self.skipTest('PowerShell unavailable; Python reuse tests still run')
        with tempfile.TemporaryDirectory(prefix='PC-4 failure 日本語 ') as temp:
            root = Path(temp)
            pc3 = root / 'PC-3'
            pc3.mkdir()
            (pc3 / 'pc3-failure.json').write_text('original PC-3 failure', encoding='utf-8')
            before = gate.source_manifest(pc3)
            pc2 = root / 'PC-2'
            for denom in (128, 256):
                for size in ('coarse', 'fine', 'finer'):
                    folder = pc2 / f'native{denom}-{size}-audit'
                    folder.mkdir(parents=True)
                    (folder / 'fold-results.json').write_text('{}', encoding='utf-8')
            # ビルド呼出しのみスタブ。数値経路へは進まず再利用検査で必ず停止する。
            bin_dir = root / 'bin'
            bin_dir.mkdir()
            cargo = bin_dir / ('cargo.cmd' if os.name == 'nt' else 'cargo')
            cargo.write_text('@exit /b 0\n' if os.name == 'nt' else '#!/bin/sh\nexit 0\n', encoding='ascii')
            cargo.chmod(0o755)
            env = dict(os.environ, PATH=str(bin_dir) + os.pathsep + os.environ['PATH'])
            output = root / 'PC-4'
            command = [shell, '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File',
                       str(EVIDENCE / 'rerun-pc.ps1'), '-Python', sys.executable,
                       '-Pc2', str(pc2), '-ResumeFrom', str(pc3), '-Output', str(output)]
            child = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, encoding='utf-8')
            self.assertNotEqual(child.returncode, 0)
            failure = gate.read_json(output / 'pc4-failure.json')
            self.assertEqual(failure['failure_stage'], 'pc3-reuse-validation', child.stdout + child.stderr)
            self.assertTrue(failure['failed_conditions'])
            self.assertIs(failure['passed'], False)
            self.assertIs(failure['nc01_complete'], False)
            self.assertIs(failure['p6_complete'], False)
            self.assertEqual(gate.source_manifest(pc3), before)
            saved = (output / 'pc4-failure.json').read_bytes()
            self.assertNotEqual(subprocess.run(command, cwd=ROOT, env=env, capture_output=True).returncode, 0)
            self.assertEqual((output / 'pc4-failure.json').read_bytes(), saved)
            self.assertFalse((output / 'independent400-coarse').exists())


if __name__ == '__main__':
    unittest.main()
