"""Read-only baseline/current probes; Windows handles are entirely fake."""
import argparse
import ast
import copy
import ctypes
import ctypes.wintypes as w
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BASE = 'd8616cefb526e70f422ec8c164afd2353c8f8149'
PAUSE = 'interim/2026-10-30-pause19-durable-SRM-mixed-elements/pause/scripts/safe-pause-latest102.py'


def source(path):
    return subprocess.run(['git', 'show', BASE+':'+path], cwd=ROOT,
                          capture_output=True, check=True).stdout


def pause_probe():
    raw = source(PAUSE)
    tree = ast.parse(raw)
    loop = next(n for n in ast.walk(tree) if isinstance(n, ast.For)
                and isinstance(n.iter, ast.Name) and n.iter.id == 'plan')
    rows = []
    for mismatch in ('image', 'creation_time'):
        for optimize in (0, 1):
            calls = []
            class Kernel:
                def OpenProcess(self, *args):
                    calls.append('OpenProcess')
                    return 1
                def QueryFullProcessImageNameW(self, handle, flags, image, size):
                    calls.append('QueryFullProcessImageNameW')
                    image.value = 'C:/unrelated.exe' if mismatch == 'image' else 'C:/owned.exe'
                    return 1
                def GetProcessTimes(self, *args):
                    calls.append('GetProcessTimes')
                    return 1  # zero FILETIME cannot match the 2026 snapshot
            namespace = dict(plan=[dict(ProcessId=100, ExecutablePath='C:/owned.exe',
                                        Created='2026-10-30T22:00:00+00:00')],
                             handles=[], kernel=Kernel(), ctypes=ctypes, w=w,
                             datetime=datetime, norm=lambda x: x.replace('\\', '/').lower())
            error = None
            try:
                exec(compile(ast.Module(body=[loop], type_ignores=[]), PAUSE, 'exec',
                             optimize=optimize), namespace)
            except AssertionError as failure:
                error = type(failure).__name__+': '+str(failure)
            rows.append(dict(mismatch=mismatch, optimize=optimize,
                             validation_block_accepted=error is None, error=error, calls=calls))
    return dict(source_sha256=hashlib.sha256(raw).hexdigest(), results=rows,
                scope='AST-extracted original held-handle validation block; no CIM, OS handle, signal, termination, or receipt operation')


def run(output):
    sys.path[:0] = [str(ROOT/'tests'), str(ROOT/'tools')]
    from test_issue_triage_followups import boundary_fixture
    from verify_soil_source_column import audit as current
    report = dict(baseline_commit=BASE, pause=pause_probe(), soil=[])
    raw = source('tools/verify_soil_source_column.py')
    report['soil_baseline_source_sha256'] = hashlib.sha256(raw).hexdigest()
    with tempfile.TemporaryDirectory() as temp:
        script = Path(temp)/'baseline_soil.py'
        script.write_bytes(raw)
        spec = importlib.util.spec_from_file_location('baseline_soil', script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        normal = boundary_fixture()
        before = module.audit(*normal)
        after = current(*normal)
        canonical = lambda value: json.dumps(value, sort_keys=True, allow_nan=False).encode()
        report['normal_report'] = dict(equal=before == after,
            baseline_sha256=hashlib.sha256(canonical(before)).hexdigest(),
            current_sha256=hashlib.sha256(canonical(after)).hexdigest())
        for defect in ('trailing_no_flow_NaN', 'negative_infinite_balance_error', 'negative_balance_error'):
            values = copy.deepcopy(normal)
            solution = values[-1]['flow_audits'][0]['solution']
            if defect == 'trailing_no_flow_NaN':
                faces = [p for p in solution['boundary_fluxes'] if p['boundary_kind'] == 'no_flow']
                faces[-1]['gradient_outflow_m3_s'] = float('nan')
            else:
                solution['balances'][0]['relative_nodal_balance_error'] = -float('inf') if 'infinite' in defect else -1.
            row = dict(defect=defect, api={})
            for label, audit in (('baseline', module.audit), ('current', current)):
                try:
                    result = audit(*values)
                    row['api'][label] = dict(accepted=True, passed=result['passed'],
                        maximum_no_flow_face_outflow_m3_s=result['frames'][0]['maximum_no_flow_face_outflow_m3_s'])
                except (AssertionError, ValueError, TypeError) as error:
                    row['api'][label] = dict(accepted=False, error=type(error).__name__+': '+str(error))
            row['cli'] = []
            for label, path in (('baseline', script), ('current', ROOT/'tools/verify_soil_source_column.py')):
                for mode in ('normal', '-O', 'env'):
                    folder = Path(temp)/(defect+'-'+label+'-'+mode)
                    folder.mkdir()
                    for name, value in zip(('mesh', 'ground', 'stages', 'stage-results'), values, strict=True):
                        (folder/(name+'.json')).write_text(json.dumps(value), encoding='utf-8')
                    destination = folder/'report.json'
                    env = dict(os.environ, PYTHONPATH=str(ROOT/'tools'))
                    env.pop('PYTHONOPTIMIZE', None)
                    if mode == 'env':
                        env['PYTHONOPTIMIZE'] = '1'
                    result = subprocess.run([sys.executable, *(['-O'] if mode == '-O' else []),
                        str(path), '--evidence', str(folder), '--output', str(destination)],
                        env=env, capture_output=True, text=True)
                    row['cli'].append(dict(version=label, mode=mode, exit_code=result.returncode,
                        report_created=destination.exists(), passed=json.loads(destination.read_text())['passed'] if destination.exists() else None,
                        stderr_last_line=result.stderr.splitlines()[-1] if result.stderr else None))
            report['soil'].append(row)
    report['scope'] = 'Analytic 12 Tet10 fixture and real checker API/CLI. Baseline checker bytes come from git show; dependencies and input fixture use the current checkout. No Rust solve, real process stop, CAD, or GUI/GPU.'
    with output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output)
