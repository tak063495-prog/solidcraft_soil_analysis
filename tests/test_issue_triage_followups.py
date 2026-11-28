"""Issue #5/#14 follow-ups: isolated entry points and external boundary data."""
import copy
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from test_p2_soil_source import fixture

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import fem_reference as fem
from joint_water_reference import Water
from verify_soil_source_column import audit

PAUSE = ROOT/'interim/2026-10-30-pause19-durable-SRM-mixed-elements/pause/scripts/safe-pause-latest102.py'


def boundary_fixture():
    """Use all twelve physical exterior faces, including ten no-flow faces."""
    values = fixture()
    mesh, ground, case, native = values
    xyz = np.asarray(mesh['nodes'])
    cells = fem.cells(mesh)
    owners = {}
    for ei in np.flatnonzero(Water(mesh, ground, case).soil):
        for face in itertools.combinations(cells[ei][:4], 3):
            owners.setdefault(tuple(sorted(map(int, face))), []).append(int(ei))
    faces = []
    for face, live in sorted(owners.items()):
        if len(live) != 1:
            continue
        ei = live[0]
        a, b, c = xyz[list(face)]
        center = (a+b+c)/3
        area = np.cross(b-a, c-a)/2
        if area @ (center-xyz[cells[ei][:4]].mean(axis=0)) < 0:
            area = -area
        kind = 'head' if np.all(xyz[list(face), 2] == 0) else 'no_flow'
        flux = float(np.array([0., 0., 2e-4*(center[2]-.5)]) @ area)
        faces.append(dict(corners=list(face), owner=ei, boundary_kind=kind,
                          gradient_outflow_m3_s=flux,
                          prescribed_outflow_m3_s=None if kind == 'head' else 0.))
    native['flow_audits'][0]['solution']['boundary_fluxes'] = faces
    return values


def cli(values, folder, mode):
    for name, value in zip(('mesh', 'ground', 'stages', 'stage-results'), values, strict=True):
        (folder/(name+'.json')).write_text(json.dumps(value), encoding='utf-8')
    output = folder/'report.json'
    env = dict(os.environ)
    env.pop('PYTHONOPTIMIZE', None)
    if mode == 'env':
        env['PYTHONOPTIMIZE'] = '1'
    command = [sys.executable, *(['-O'] if mode == '-O' else []),
               str(ROOT/'tools/verify_soil_source_column.py'),
               '--evidence', str(folder), '--output', str(output)]
    result = subprocess.run(command, env=env, capture_output=True, text=True)
    return result, output


class SoilBoundaryFlow(unittest.TestCase):
    def test_normal_complete_boundary_inventory_in_all_modes(self):
        values = boundary_fixture()
        faces = values[-1]['flow_audits'][0]['solution']['boundary_fluxes']
        self.assertEqual(len(faces), 12)
        self.assertEqual(sum(p['boundary_kind'] == 'no_flow' for p in faces), 10)
        expected = audit(*values)
        self.assertTrue(expected['passed'])
        self.assertEqual(expected['frames'][0]['maximum_no_flow_face_outflow_m3_s'], 0.)
        for mode in ('normal', '-O', 'env'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as folder:
                result, output = cli(values, Path(folder), mode)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(json.loads(output.read_text()), expected)

    def test_first_and_last_nonfinite_non_scalar_or_mismatched_flow_are_rejected(self):
        base = boundary_fixture()
        for kind in ('head', 'no_flow'):
            for end in (0, -1):
                for bad in (float('nan'), float('inf'), -float('inf'), [], [0.], [[0.]], 1.):
                    values = copy.deepcopy(base)
                    faces = values[-1]['flow_audits'][0]['solution']['boundary_fluxes']
                    selected = [p for p in faces if p['boundary_kind'] == kind][end]
                    selected['gradient_outflow_m3_s'] = bad
                    with self.subTest(kind=kind, end=end, bad=bad):
                        with self.assertRaises((ValueError, AssertionError, TypeError)):
                            audit(*values)
                    for mode in ('normal', '-O', 'env'):
                        with self.subTest(kind=kind, end=end, bad=bad, mode=mode), tempfile.TemporaryDirectory() as folder:
                            result, output = cli(values, Path(folder), mode)
                            self.assertNotEqual(result.returncode, 0, result.stdout)
                            self.assertFalse(output.exists(), result.stdout)

    def test_reported_balance_error_must_be_finite_scalar_and_nonnegative(self):
        for bad in (float('nan'), float('inf'), -float('inf'), -1., [], [0.], None):
            values = boundary_fixture()
            values[-1]['flow_audits'][0]['solution']['balances'][0]['relative_nodal_balance_error'] = bad
            with self.subTest(bad=bad):
                with self.assertRaises((ValueError, AssertionError, TypeError)):
                    audit(*values)
            for mode in ('normal', '-O', 'env'):
                with self.subTest(bad=bad, mode=mode), tempfile.TemporaryDirectory() as folder:
                    result, output = cli(values, Path(folder), mode)
                    self.assertNotEqual(result.returncode, 0, result.stdout)
                    self.assertFalse(output.exists(), result.stdout)


class SafePauseEntry(unittest.TestCase):
    def test_normal_help_remains_read_only(self):
        env = dict(os.environ)
        env.pop('PYTHONOPTIMIZE', None)
        result = subprocess.run([sys.executable, str(PAUSE), '--help'], env=env,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--execute', result.stdout)

    def test_optimized_execute_is_refused_before_process_or_filesystem_work(self):
        # Even a regressed guard cannot reach CIM, native handles, or a real job.
        harness = """import ctypes,runpy,subprocess,sys
from pathlib import Path
def forbidden(*a, **kw):
    raise RuntimeError('unexpected process or filesystem operation')
subprocess.run = subprocess.check_output = subprocess.Popen = forbidden
ctypes.WinDLL = forbidden
Path.mkdir = Path.touch = Path.write_text = Path.write_bytes = forbidden
script = sys.argv[1]
sys.argv = [script, '--execute', '--usage-percent', '81']
runpy.run_path(script, run_name='__main__')
"""
        for mode in ('-O', 'env'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as folder:
                staged = Path(folder)/'root/pause/scripts/safe-pause-latest102.py'
                staged.parent.mkdir(parents=True)
                staged.write_bytes(PAUSE.read_bytes())
                before = sorted(str(p.relative_to(folder)) for p in Path(folder).rglob('*'))
                env = dict(os.environ)
                env.pop('PYTHONOPTIMIZE', None)
                if mode == 'env':
                    env['PYTHONOPTIMIZE'] = '1'
                result = subprocess.run([sys.executable, *(['-O'] if mode == '-O' else []),
                                         '-c', harness, str(staged)], cwd=folder, env=env,
                                        capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('safe-pause requires Python optimization disabled', result.stderr)
                self.assertNotIn('unexpected process or filesystem operation', result.stderr)
                self.assertEqual(result.stdout, '')
                self.assertEqual(sorted(str(p.relative_to(folder)) for p in Path(folder).rglob('*')), before)


if __name__ == '__main__':
    unittest.main()
