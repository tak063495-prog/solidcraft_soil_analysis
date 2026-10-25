"""Create explicit domain-change examples; never overwrite existing inputs."""
import json
from pathlib import Path

def examples(root):
    for model in ('mc','mcc'):
        base=json.loads((root/f'examples/one-way-{model}.json').read_text(encoding='utf-8'))
        base['flow']['analysis']['periods']=[dict(name='before',end_time_s=200.,boundaries=[
            dict(face_set='bottom',condition=dict(kind='head',value=dict(head_m=0.))),
            dict(face_set='top',condition=dict(kind='atmospheric_no_ponding',rainfall_vertical_m_s=0.))]),
            dict(name='after',end_time_s=400.,boundaries=[
            dict(face_set='bottom',condition=dict(kind='head',value=dict(head_m=0.))),
            dict(face_set='floor',condition=dict(kind='atmospheric_no_ponding',rainfall_vertical_m_s=0.))])]
        base['exterior_periods']=[dict(name=name,reservoirs=[]) for name in ('before','after')]
        case=dict(schema='soil-fem-one-way-construction/1',one_way=base,excavations=[dict(name='cut',before_period='after',bodies=['UpperClay'],water_boundary_policy='explicit_period_boundaries',cut_traction_policy='explicit_period_reservoirs')])
        with (root/f'examples/one-way-construction-{model}.json').open('x',encoding='utf-8') as stream:json.dump(case,stream,indent=2)

if __name__=='__main__':examples(Path(__file__).resolve().parents[1])
