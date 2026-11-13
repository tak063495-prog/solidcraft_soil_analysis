"""Read-only extraction of original five failed gates; never edits evidence."""
import argparse
import hashlib
import json
from pathlib import Path

BASE = '702a19a7ce189f1173caed727ee2a0059298a607'
OLD = 'interim/2026-10-29-native-kernel-and-rain-gates/rain-original/'
PATHS = [
    'interim/2026-11-16-issue-triage/todo.ja.md', 'docs/P6_STATUS.ja.md',
    'docs/RICHARDS.ja.md', 'docs/ONE_WAY.ja.md',
    OLD+'scope.json', OLD+'reference-case.json', OLD+'comparison.json', OLD+'compare.py',
    'interim/2026-10-29-epoch-export-rain-accuracy/rain-refined-comparison.json',
    'interim/2026-10-30-rain-refinement-and-finer14909/README.ja.md',
    'interim/2026-10-30-rain-refinement-and-finer14909/rain/comparison.json']


def extract(root):
    root = Path(root)
    read = lambda p: json.loads((root/p).read_text(encoding='utf-8-sig'))
    scope, case, comparison = [read(OLD+p) for p in ('scope.json', 'reference-case.json', 'comparison.json')]
    physical = dict(geometry=scope['geometry'], rain_patch_m=scope['rain_patch_m'],
                    rain_area_m2=scope['rain_area_m2'], materials=case['materials'],
                    initial_total_head_m=case['initial_head']['head_m'],
                    initial_pressure_head='psi(z)=-z; z=0 bottom, z=1 top',
                    initial_saturation='Sr(z)=theta_VG(-z)/theta_s (not effective saturation)',
                    boundary_conditions='bottom total H=0; sides no-flow; top atmospheric_no_ponding; patch rain, dry remainder 0',
                    physical_periods=scope['physical_periods'],
                    native_subdivisions=scope['native_subdivisions'],
                    native_tet4_counts=[384, 3072, 24576],
                    native_dt_s=scope['native_time_step_s'],
                    reference_subdivisions=scope['independent_subdivisions'],
                    reference_dt_s=scope['independent_time_step_s'],
                    mechanics=False, surface_storage_or_runoff_dynamics=False)
    rows = []
    for index, check in enumerate(comparison['checks']):
        if check['passed']:
            continue
        stage = check['name'].split()[0]
        time_s = {'rain': 3600., 'pause': 7200., 'rain_again': 10800.}[stage]
        ratio = 'ratio' in check['name']
        metric = ('RMS difference(8,16)/RMS difference(4,8), RT0 six-tet block averages trilinearly sampled at common FV32 centres'
                  if ratio else 'native16 cube-mean head vs trilinear FV32 at native16 centres'
                  if 'head' in check['name'] else 'absolute native16 minus FV32 total storage offset')
        row = dict(id=f'OLD-RAIN-{len(rows)+1:02}', name=check['name'], physical=physical,
                   time_s=time_s, comparison_quantity=metric,
                   unit='1' if ratio else 'm' if 'head' in check['name'] else 'm3',
                   original_criterion='finite(actual) and actual < upper_bound',
                   original_upper_bound=check['upper_bound'], original_recorded_actual=check['actual'],
                   original_passed=check['passed'], current_status='未達（旧記録を保持）',
                   evidence_pointer=f'/checks/{index}',
                   failure_reason='元の厳密な上限を超える。離散化の根本原因はこの記録だけでは未確認。',
                   independent_same_condition_value=None, accepted=None,
                   independent_reason='未定義：有限正方形の降雨パッチと二横方向流を同時に保存する1D/2D縮約なし。新規3D求解は未実行。',
                   missing_input_values=[],
                   not_confirmed=['元3Dと同じ新規独立参照値', '不合格の数値アルゴリズム上の根本原因'])
        rows.append(row)
    if len(rows) != 5 or comparison['passed'] is not False:
        raise ValueError('original five failure record differs')
    fingerprints = {p: {'sha256': hashlib.sha256((root/p).read_bytes()).hexdigest(),
                        'bytes': (root/p).stat().st_size, 'commit': BASE} for p in PATHS}
    later = []
    for p in PATHS:
        if 'comparison.json' in p and p != OLD+'comparison.json':
            data = read(p)
            later.append({'path': p, 'passed': data.get('passed'),
                          'check_count': len(data.get('checks', [])),
                          'meaning': '後続細分化の別範囲。旧5不合格の置換禁止'})
            if p == 'interim/2026-10-30-rain-refinement-and-finer14909/rain/comparison.json':
                by_name = {c['name']: c for c in data['checks']}
                for row in rows:
                    row['later_bounded_record'] = dict(
                        path=p, check=by_name.get(row['name']),
                        scope='native 8/16/24, finest 82944 Tet4 vs FV64; different refinement window',
                        replaces_original_failure=False)
    return dict(base_commit=BASE, old_rows=rows, original_predeclared_gates=scope['predeclared_gates'],
                later_bounded_records=later, source_fingerprints=fingerprints,
                p6_03_complete=False, old_failure_count=5)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    Path(args.out).write_text(json.dumps(extract(args.root), ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
