"""Read historical evidence and classify its relationship to new references.

Never imports/calls any application/reference solver and never edits evidence.
"""
import argparse
import hashlib
import json
from pathlib import Path


EVIDENCE = [
    'docs/P6_STATUS.ja.md', 'docs/STRUCTURE_WATER_JOINT_REMAINDERS.ja.md',
    'docs/SOLID_PILE.ja.md', 'interim/2026-11-16-issue-triage/todo.ja.md',
    'examples/solid-pile.json','examples/ground-solid-pile.json',
    'examples/stages-solid-pile.json','examples/stages-solid-pile-body-planes.json',
    'interim/2026-10-26-body-pile/reports/solid-pile-first-pair.json',
    'interim/2026-10-26-body-pile/reports/solid-pile-refined-pair.json',
    'tools/verify_solid_pile.py', 'scripts/Invoke-SolidPile.ps1',
]


def diagnose(root, small_3d=None):
    root=Path(root)
    fingerprints={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in EVIDENCE}
    binding=json.loads(Path(__file__).with_name('source-bindings.json').read_text(encoding='utf-8'))
    if fingerprints != binding['files_sha256']:
        raise ValueError('historical source bytes differ from pinned source bindings; no implicit comparison')
    first=json.loads((root/EVIDENCE[8]).read_text(encoding='utf-8-sig'))
    refined=json.loads((root/EVIDENCE[9]).read_text(encoding='utf-8-sig'))
    values=refined['three_mesh_loading_settlement_mm']
    if values[:2] != first['head_area_mean_loading_settlement_mm']:
        raise ValueError('historical records disagree')
    if first['spatial_convergence_accepted'] is not False or refined['spatial_convergence_accepted'] is not False:
        raise ValueError('historical judgement changed; review manually')
    relative=[abs(b-a)/abs(b) for a,b in zip(values[:-1],values[1:])]
    historical=[first['coarse_to_fine_settlement_relative_difference'],refined['current_pair_settlement_relative_difference']]
    if any(abs(a-b)>1e-14 for a,b in zip(relative,historical)):
        raise ValueError('historical metric definition does not match record')
    rows=[]
    if small_3d is not None:
        if small_3d['schema']!='pile-independent-small-3d/1' or small_3d['accepted'] is not None:
            raise ValueError('require nonaccepted standalone 3D report')
        computed=[r['head_additional_settlement_mm'] for r in small_3d['series']]
        for count,value in zip(refined['three_mesh_tet10_counts'],values):
            rows.append(dict(historical_tet10=count,historical_settlement_mm=value,
                             sampled_hex8_value_range_position='含まれる' if min(computed)<=value<=max(computed) else '外れる',
                             difference_from_finest_hex8_mm=value-computed[-1],
                             independent_reference_judgement='判定不能',accepted=None))
    return dict(schema='pile-historical-diagnosis/1',source_commit='7ae0a2a669e0e1df93f3dd9791bbeda056334c3e',
                historical_source_sha256=fingerprints,
                historical_tet10_counts=refined['three_mesh_tet10_counts'],
                historical_settlement_mm=values, historical_relative_difference=historical,
                recomputed_relative_difference=relative,
                metric='abs(finer head area-mean load-increment settlement - coarser)/abs(finer)',
                metric_is_exact_solution_error=False,spatial_threshold='未定義',accepted=None,
                historical_spatial_convergence_accepted=False,
                axial_reference_position='判定不能: shaft/tip stiffness and square finite-box mapping 未確認',
                winkler_reference_position='判定不能: lateral spring and finite-box mapping 未確認',
                same_physics_small_hex8=rows,
                sampled_range_is_error_bound=False,
                missing=['historical all-node/GP output in public first/refined summary records',
                         'independent error-bounded/converged continuum solution',
                         'verified axial/lateral spring mapping and tip stiffness',
                         'pile spatial-convergence tolerance definition'],
                error_sources=['historical Tet10 vs own Hex8 spatial discretization and nonuniform mesh shape',
                               'unresolved interface/corner gradients and finite side/bottom boundaries',
                               '1D springs omit continuum interaction, finite section and Poisson/shear effects'],
                p6_complete=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    parser.add_argument('--small-3d',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=None if args.small_3d is None else json.loads(args.small_3d.read_text(encoding='utf-8'))
    result=diagnose(args.root,report)
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False);stream.write('\n')


if __name__=='__main__':main()
