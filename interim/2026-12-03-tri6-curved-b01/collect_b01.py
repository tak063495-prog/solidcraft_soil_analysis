"""Summarize solved paths; retain false/null and undefined gates."""
import argparse
import copy
from pathlib import Path
from b01_cases import read, write


def run(work, out):
    fresh = work/'b01-fresh'
    maximum, cases = {}, []
    for original in read(fresh/'measurements.json')['cases']:
        row = copy.deepcopy(original)
        result = read(fresh/row['case']/'native/run.json')
        attempts = result['attempts']+result['geostatic']['attempts']
        accepted = [a for a in attempts if a['accepted']]
        linear = [s['true_relative_residual'] for a in accepted for s in a['linear_solves']]
        row['solver'] = dict(complete=result['complete'], accepted=len(accepted),
            rejected=sum(not a['accepted'] for a in attempts),
            maximum_accepted_force_residual=max(a['relative_force_residual'] for a in accepted),
            maximum_true_linear_residual=max(linear, default=0.))
        cmp = row['comparison']
        row['existing_limited_gates_passed'] = (cmp['existing_field_gates_passed']
            and cmp['maximum_native_reassembled_residual'] < 1e-10
            and cmp['maximum_independent_reassembled_residual'] < 1e-10
            and row['solver']['maximum_accepted_force_residual'] < 1e-8
            and row['solver']['maximum_true_linear_residual'] < 1e-10)
        for key, value in cmp['maximum_by_field'].items():
            item = maximum.setdefault(key, dict(maximum_absolute=0., relative_l2=0., components_compared=0,
                                                gate=value['gate'], passed=value['passed']))
            for metric in ('maximum_absolute', 'relative_l2'):
                item[metric] = max(item[metric], value[metric])
            item['components_compared'] += value['components_compared']
            if value['passed'] is False:
                item['passed'] = False
        cases.append(row)
    write(out/'measurements.json', dict(passed=False, b01_complete=False, all_field_numerical_acceptance=None,
        cases=cases, maximum_by_field=maximum,
        existing_limited_gates_passed=all(r['existing_limited_gates_passed'] for r in cases),
        accepted_point_pairs=sum(r['comparison']['accepted_points'] for r in cases),
        maximum_native_reassembled_residual=max(r['comparison']['maximum_native_reassembled_residual'] for r in cases),
        maximum_independent_reassembled_residual=max(r['comparison']['maximum_independent_reassembled_residual'] for r in cases),
        maximum_true_linear_residual=max(r['solver']['maximum_true_linear_residual'] for r in cases),
        maximum_accepted_force_residual=max(r['solver']['maximum_accepted_force_residual'] for r in cases)))
    for source, name in [('b01-probes/scope.json', 'scope-refusals.json'),
                         ('b01-mutations/mutations.json', 'checkpoint-mutations.json'),
                         ('b01-fresh/restarts.json', 'restart-audit.json'),
                         ('b01-curve/curved-comparison.json', 'curved-comparison.json'),
                         ('b01-cli-controls/cli-controls.json', 'production-cli-controls.json')]:
        write(out/name, read(work/source))
    decisions = [
        ('平面・直線Tri6のe=0多面体、両端地盤付着の円筒2/4/8分割', '限定受入の再確認：今回の乾燥弾性・直接入力・小型範囲のみ'),
        ('平面・直線Tri6のe=0古典屋根2/4分割', '限定受入の再確認：実装/全場一致のみ。曲面h収束の数値受入は未達'),
        ('面内で曲がった中間節点／曲線辺', '未達／未対応（現行拒否を保持）'),
        ('面外中間節点・要素内で法線が変わる二次曲面・真曲面パッチ', '未達／未対応（現行拒否を保持）'),
        ('真曲面＋正負オフセット', '未達／未対応（現行拒否を保持）'),
        ('共有IDの非ゼロ腕で法線/e不一致・範囲外e・宣言偽装', '未達／未対応（現行拒否を保持）'),
        ('平面offsetの既存限定範囲、旧mindlin/cst/既定のeなし', '既存動作を保持。今回の4 CLI回帰対照も完走'),
        ('曲面Tri6仕様/適合製造解/曲面MITC/drilling/円筒・半球の本体h系列', '未達／未対応（現行拒否を保持）。古典屋根の多面体比較で代替しない'),
        ('コンパクト/Native曲面準備・GUI/GPU・実CAD', '未達／未対応。既存拒否を保持、PC手順のみ'),
    ]
    write(out/'acceptance.json', dict(passed=False, b01_complete=False, curved_shell_acceptance=None,
        combinations=[dict(combination=k, judgment=v) for k, v in decisions],
        no_production_scope_or_default_change=True, all_other_remainders_preserved=True,
        old_B02_B03_failures_and_NC_records_preserved=True))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('work', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    run(a.work, a.output)
