"""PC rate-arc の限定ゲート。数値求解を行わず、診断の passed=false を保つ。"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

PC3_HEAD = '8342762bd714e33523ce5fef772572699fc21881'


def check_rate_arc(report, samples):
    """PowerShell の配列列挙・暗黙型変換に依存しない、固定128点の検査。"""
    failures = []
    if type(report) is not dict:
        return ['report: JSON object required']
    required = ('equilibrium_samples', 'stop_reason',
                'maximum_bordered_true_linear_residual', 'passed')
    failures.extend(f'report.{key}: missing' for key in required if key not in report)
    count = report.get('equilibrium_samples')
    if type(count) is not int or count != 128:
        failures.append('report.equilibrium_samples: integer 128 required')
    if type(samples) is not list or len(samples) != 128:
        failures.append('samples: top-level array of exactly 128 objects required')
    elif any(type(row) is not dict for row in samples):
        failures.append('samples: each row must be an object (nested arrays rejected)')
    if 'stop_reason' in report and report['stop_reason'] is not None:
        failures.append('report.stop_reason: JSON null required')
    if report.get('passed') is not False:
        failures.append('report.passed: boolean false required')
    residual = report.get('maximum_bordered_true_linear_residual')
    if (type(residual) not in (int, float) or not 0 <= residual <= 1e-10
            or not math.isfinite(residual)):
        failures.append('report.maximum_bordered_true_linear_residual: finite number in [0, 1e-10] required')
    return failures


def read_json(path):
    # Windows PowerShell の Out-File UTF8 は BOM を付ける。
    def number(value):
        parsed = float(value)
        if not math.isfinite(parsed):
            raise ValueError('non-finite JSON number: ' + value)
        return parsed
    def invalid(value):
        raise ValueError('invalid JSON constant: ' + value)
    return json.loads(Path(path).read_text(encoding='utf-8-sig'),
                      parse_float=number, parse_constant=invalid)


def inspect_arc(arc):
    failures = []
    values = {}
    for name in ('report', 'samples'):
        try:
            values[name] = read_json(arc / f'{name}.json')
        except (OSError, ValueError) as exc:
            failures.append(f'{name}.json: {exc}')
    if not failures:
        failures = check_rate_arc(values['report'], values['samples'])
    report = values.get('report')
    return dict(accepted=not failures, passed=False, nc01_complete=False,
                p6_complete=False, gate='rate_arc_128', failures=failures,
                arc=str(arc.resolve()), report=report,
                samples_type=type(values.get('samples')).__name__,
                samples_count=(len(values['samples']) if type(values.get('samples')) is list else None))


def source_manifest(root):
    """再利用元の全ファイルを読み取り専用で指紋化する。失敗記録も含む。"""
    paths = [root / 'head.txt', root / 'reference400-command.json', root / 'rate-arc400-command.json']
    paths += sorted(p for folder in ('reference400', 'rate-arc400')
                    for p in (root / folder).rglob('*') if p.is_file())
    paths += sorted(root.glob('*failure*'))
    return {str(p.relative_to(root)): file_sha256(p)
            for p in paths if p.is_file()}


def file_sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def inspect_reuse(root, repo):
    """PC-3 の reference400/rate-arc400 だけを検証。別HEADのcheckpointは再開しない。"""
    result = inspect_arc(root / 'rate-arc400')
    result['gate'] = 'pc3_reuse'
    failures = result['failures']
    try:
        head = (root / 'head.txt').read_text(encoding='utf-8-sig').strip()
        resolved = subprocess.check_output(['git', 'rev-parse', '--verify', head + '^{commit}'],
                                           cwd=repo, text=True).strip()
        if resolved != head or head != PC3_HEAD:
            raise ValueError('PC-3 head: expected ' + PC3_HEAD)
        result['source_head'] = head
        for name, script, options in (
            ('reference400', 'tools/low_preload_nc_reference.py', {'--steps': '400', '--method': 'rk4'}),
            ('rate-arc400', 'interim/2026-11-09-nc-pseudo-arclength/independent_arc.py',
             {'--count': '128', '--step-size': '0.00002'}),
        ):
            record = read_json(root / f'{name}-command.json')
            if type(record.get('return_code')) is not int or record['return_code'] != 0 or record.get('head') != head:
                failures.append(f'{name}-command: rc=0 and source head required')
            args = record.get('arguments')
            if type(args) is not list or not args or args[0].replace('\\', '/') != script:
                failures.append(f'{name}-command: unexpected script/arguments')
                continue
            for key, value in options.items():
                if args.count(key) != 1 or args.index(key) + 1 >= len(args) or args[args.index(key) + 1] != value:
                    failures.append(f'{name}-command: expected {key} {value}')
        reference = read_json(root / 'reference400/reference.json')
        if (reference.get('last_accepted', {}).get('name') != 'excavate_upper'
                or reference.get('steps') != 400 or reference.get('integration_method') != 'rk4'
                or reference.get('passed') is not False or reference.get('input_only') is not True):
            failures.append('reference400: excavate_upper / steps=400 / rk4 / passed=false / input_only=true required')
        for name in ('mesh', 'ground', 'stages', 'unload'):
            if type(read_json(root / f'reference400/{name}.json')) is not dict:
                failures.append(f'reference400/{name}.json: object required')
        # 既存生成器の入力と全サンプルの形状を検査するだけで、再求解しない。
        import numpy as np
        from probe_low_preload_nc import expected_inputs
        mesh, ground, case = expected_inputs()
        for name, expected in zip(('mesh', 'ground', 'stages'), (mesh, ground, case)):
            if read_json(root / f'reference400/{name}.json') != expected:
                failures.append(f'reference400/{name}.json: differs from unchanged input generator')
        samples = read_json(root / 'rate-arc400/samples.json')
        active = reference['last_accepted']['active_elements']
        if len(active) != len(mesh['tetrahedra']) or any(type(v) is not bool for v in active):
            failures.append('reference400.last_accepted.active_elements: invalid mask')
        for index, row in enumerate(samples):
            if row.get('step') != index or row.get('input_only') is not True or row.get('step_size_m') != 2e-5:
                failures.append(f'samples[{index}]: step/input_only/step_size_m mismatch')
            for key, shape in (('displacement_m', (len(mesh['nodes']), 3)),
                               ('material_state', (len(mesh['tetrahedra']), 27, 14)),
                               ('tangent', (55,))):
                array = np.asarray(row[key], dtype=float)
                if array.shape != shape or not np.isfinite(array).all():
                    failures.append(f'samples[{index}].{key}: shape or non-finite value')
            if type(row.get('progress')) not in (int, float) or not math.isfinite(row['progress']):
                failures.append(f'samples[{index}].progress: finite number required')
            if type(row.get('contacts')) is not dict:
                failures.append(f'samples[{index}].contacts: object required')
        if result['report'].get('diagnostic_status') != 'bounded_samples':
            failures.append('rate-arc400.report.diagnostic_status: bounded_samples required')
        controls = read_json(root / 'rate-arc400/controls.json')
        for key, expected in (('count', 128), ('step_size_m', 2e-5), ('input_only', True),
                              ('passed', False), ('existing_true_linear_tolerance', 1e-10)):
            if key not in controls or type(controls[key]) is not type(expected) or controls[key] != expected:
                failures.append(f'rate-arc400.controls.{key}: expected {expected!r}')
        # 前史・rateの計算ソース/入力が旧HEADと同一であることを保守的に要求。
        # 新しい検査器だけを除外する。無関係なPython変更でも停止側に倒す。
        tracked = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', head], cwd=repo, text=True).splitlines()
        paths = [p for p in tracked if
                 (p.startswith('tools/') and p.endswith('.py')) or
                 p.startswith('interim/2026-11-09-nc-pseudo-arclength/') and p.endswith('.py') or
                 p.startswith('examples/') and p.endswith('.json')]
        changed = []
        for path in paths:
            original = subprocess.check_output(['git', 'show', f'{head}:{path}'], cwd=repo)
            current = repo / path
            if not current.is_file() or current.read_bytes() != original:
                changed.append(path)
        if changed:
            failures.append('source/input differs from PC-3 head; reuse refused: ' + ', '.join(changed))
        result['source_files_checked'] = len(paths)
        result['source_manifest'] = source_manifest(root)
        result['read_only_source'] = str(root.resolve())
    except (OSError, ValueError, KeyError, TypeError, AttributeError, ImportError, subprocess.CalledProcessError) as exc:
        failures.append(f'PC-3 reuse validation: {exc}')
    result['accepted'] = not failures
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--arc', type=Path)
    source.add_argument('--resume-from', type=Path)
    source.add_argument('--verify-reuse', type=Path)
    parser.add_argument('--repo', type=Path, default=Path.cwd())
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.verify_reuse:
        saved = read_json(args.verify_reuse)
        root = Path(saved['read_only_source'])
        failures = [] if source_manifest(root) == saved['source_manifest'] else ['PC-3 source changed after validation']
        result = dict(accepted=not failures, failures=failures, passed=False, nc01_complete=False, p6_complete=False)
    else:
        result = inspect_reuse(args.resume_from.resolve(), args.repo.resolve()) if args.resume_from else inspect_arc(args.arc)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, ensure_ascii=False, allow_nan=False, indent=2)
        stream.write('\n')
    print(json.dumps(result, ensure_ascii=True, allow_nan=False))
    return 0 if result['accepted'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
