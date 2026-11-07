"""Index frozen MCC comparison reports without solving or re-integrating fields.

This evidence inventory remains separate from the dry-MC ratio-two series.
Missing or modified records are errors; they never trigger a solver fallback.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORDS = (
    ('changing_wall',
     'interim/2026-10-27-compact-pit/reports/coupled274-945-spatial-not-accepted.json',
     '95a88691c4dc156864f49c0997c3fc594bbc7d778a9d4c20ba4094a5fc88303d'),
    ('changing_wall',
     'interim/2026-10-27-compact-pit/reports/coupled945-1324-spatial-not-accepted.json',
     'cd4343ca8d0c60dce21984cc340302d02e483f17c2ea8014a428e90c61c92f31'),
    ('changing_wall',
     'interim/2026-10-28-mcc-verification/reports/uniform1324-2262-overlap.json',
     '1f192379869bd065d6549b5dc895f424e976dc9251d22f69111ade30417bba2e'),
    ('fixed_wall',
     'interim/2026-10-28-completed-studies/reports/fixed-wall2262-2865-comparison.json',
     'fbd84db91ee43e7b1586fa8206fd55449390e4a13cb571affe627832009649d6'),
    ('fixed_wall',
     'interim/2026-10-29-fixed-wall-5144/reports/comparison.json',
     '6b485f6aa0a449283f2926036e93ba42f4724574caebf2eb79459f95308964fe'),
    ('fixed_wall',
     'interim/2026-10-29-fixed-wall8061-native-epoch/comparison-corrected-labels.json',
     'd544e64995236a9dbf010e73022349bd33f5bde047d5bb05337e475d86c0130d'),
    ('fixed_wall',
     'interim/2026-10-30-offline-setup-and-rain-restarts/spatial/original8061-to14909-full-comparison.json',
     '3d9f240e232dcbb434c4b6cc5f16d4c157e9b807bceec483ca4dae13929747b3'),
)


def collect(root=ROOT):
    rows=[]
    for group,name,digest in RECORDS:
        path=Path(root)/name
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'frozen MCC comparison changed: {name}')
        source=json.loads(path.read_text(encoding='utf-8-sig'))
        if (not source['comparison_recorded'] or source['p6_complete']
                or source.get('spatial_convergence_accepted',False)
                or source.get('spatial_convergence_verified',False)):
            raise ValueError(f'expected recorded, unaccepted MCC comparison: {name}')
        rows.append(dict(group=group,source=name,source_sha256=digest,
            left_elements=source['left_elements'],right_elements=source['right_elements'],
            wall_fixed=source.get('wall_spatial_topology_and_sections_and_loads_preserved',False),
            frames=[dict(stage=f.get('name',f.get('stage')),
                         from_seed=f['from_seed'],from_construction_origin=f['from_construction_origin'],
                         stage_increment=f['stage_increment']) for f in source['frames']],
            spatial_convergence_accepted=False))
    last=[row['frames'][-1]['stage_increment'] for row in rows[-2:]]
    return dict(schema='soil-fem-mcc-spatial-record-inventory/1',
        source_only=True,new_solves=0,new_field_integrations=0,rows=rows,
        fixed_wall_5144_8061_14909_increment_difference_non_decreasing=(
            last[1]['relative_l2'] >= last[0]['relative_l2']
            and last[1]['rms_difference_m'] >= last[0]['rms_difference_m']),
        p6_complete=False,general_3d_spatial_convergence_demonstrated=False,
        order=None,estimate=None,
        limitations=['Mesh counts are not a physical h parameter; no Richardson fit from element counts.',
                     'Changing-wall, fixed-wall and dry-MC in-plane/depth studies remain separate.',
                     'Frozen published reports are indexed, not a new audit of unavailable raw multi-GB results.',
                     'Original stale labels and their correction record are preserved.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    report=collect()
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(report,stream,indent=2,allow_nan=False)


if __name__ == '__main__':
    main()
