"""Summarize real solver output without inventing reaction/work/convergence gates."""
import argparse,copy
from pathlib import Path
from b03_cases import read,write


def run(formal,scope,supplemental,mutations,out):
    out.mkdir(parents=True,exist_ok=True);measure=read(formal/'measurements.json');maximum={};cases=[]
    for row in measure['cases']:
        row=copy.deepcopy(row);f=formal/row['case'];r=read(f/'native/run.json')
        attempts=r['attempts']+r['geostatic']['attempts'];accept=[a for a in attempts if a['accepted']]
        linear=[s['true_relative_residual'] for a in accept for s in a['linear_solves']]
        solver=dict(complete=r['complete'],accepted=len(accept),rejected=sum(not a['accepted'] for a in attempts),
             maximum_accepted_force_residual=max(a['relative_force_residual'] for a in accept),maximum_true_linear_residual=max(linear,default=0.))
        cmp=row['comparison'];gates=cmp['existing_field_gates_passed'] and cmp['maximum_native_reassembled_residual']<1e-10 and cmp['maximum_independent_reassembled_residual']<1e-10 and solver['maximum_accepted_force_residual']<1e-8 and solver['maximum_true_linear_residual']<1e-10
        row.update(solver=solver,existing_limited_gates_passed=gates,all_field_acceptance=None,b03_complete=False)
        for k,v in cmp['maximum_by_field'].items():
            item=maximum.setdefault(k,dict(maximum_absolute=0.,relative_l2=0.,components_compared=0,gate=v['gate'],passed=v['passed']))
            item['maximum_absolute']=max(item['maximum_absolute'],v['maximum_absolute']);item['relative_l2']=max(item['relative_l2'],v['relative_l2']);item['components_compared']+=v['components_compared']
            if v['passed'] is False:item['passed']=False
        cases.append(row)
    groups={}
    for kind in sorted({c['kind'] for c in cases}):
        rows=[c for c in cases if c['kind']==kind];groups[kind]=dict(cases=[r['case'] for r in rows],accepted_points=sum(r['comparison']['accepted_points'] for r in rows),
             existing_limited_gates_passed=all(r['existing_limited_gates_passed'] for r in rows))
    root=dict(passed=False,b03_complete=False,all_field_numerical_acceptance=None,cases=cases,families=groups,
              maximum_by_field=maximum,accepted_point_pairs=sum(c['comparison']['accepted_points'] for c in cases),
              maximum_native_reassembled_residual=max(c['comparison']['maximum_native_reassembled_residual'] for c in cases),
              maximum_independent_reassembled_residual=max(c['comparison']['maximum_independent_reassembled_residual'] for c in cases),
              maximum_true_linear_residual=max(c['solver']['maximum_true_linear_residual'] for c in cases),
              maximum_accepted_force_residual=max(c['solver']['maximum_accepted_force_residual'] for c in cases))
    write(out/'measurements.json',root);write(out/'restart-audit.json',read(formal/'restarts.json'))
    write(out/'scope-refusals.json',read(scope));write(out/'supplemental-refusals.json',read(supplemental));write(out/'checkpoint-mutations.json',read(mutations))
    decision=[
        dict(combination='shared-ID finite-angle nonzero-arm fold, one-zero-arm fold',judgment='未達／未対応（現行拒否を保持）'),
        dict(combination='shared-ID equal-normal e step, only normal reversed, only e reversed',judgment='未達／未対応（現行拒否を保持）'),
        dict(combination='zero-arm e=0 90-degree two-plane patch',judgment='限定受入の再確認：直線平面要素・乾燥弾性MC・この小型直接入力のみ',family='kink-zero'),
        dict(combination='coplanar simultaneous normal/e reversal with identical e*n',judgment='限定受入の再確認：同じ剛体腕を持つ平面パッチのみ',family='normal-sign-pair'),
        dict(combination='different structural IDs on same soil node; independent rotations; e step / nonzero-arm fold',judgment='独立回転の明示モデルだけの限定確認。連続段差/折れ目のB-03は未達',families=['step-split','kink-split']),
        dict(combination='two parallel planes z=.5/1 without shared structural nodes',judgment='独立平板群だけの限定確認。段差継手の受入は未達',family='height-step'),
        dict(combination='nonplanar/curved midside node or attachment mismatch',judgment='未達／未対応（現行拒否を保持）'),
        dict(combination='compact offset surfaces / Native preparation and GUI/GPU / real CAD',judgment='未達／未対応（現行拒否を保持）。GUI/GPU/実CADはPC手順のみ'),
    ]
    write(out/'acceptance.json',dict(passed=False,b03_complete=False,all_field_numerical_acceptance=None,combinations=decision,family_measurements=groups,
          no_production_scope_or_default_change=True,all_other_remainders_preserved=True,
          inherited_B02_768_Tet10_reassembly_and_negative_spatial_orders_preserved=True))
if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('formal','scope','supplemental','mutations','output'):p.add_argument(k,type=Path)
    a=p.parse_args();run(a.formal,a.scope,a.supplemental,a.mutations,a.output)
