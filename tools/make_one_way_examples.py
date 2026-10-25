"""Deterministic physical-time soil/water examples; all output is explicit UTF-8."""
import copy
import json
from pathlib import Path
root=Path(__file__).resolve().parents[1]


def build(mcc=False):
    ground=json.loads((root/('examples/ground-seepage-column-mcc.json' if mcc else 'examples/ground-water-column-mc.json')).read_text(encoding='utf-8'))
    ground['initialization']['water_level_m']=None
    for m in ground['materials']:m['unit_weight_kn_m3']=16.
    mechanics=json.loads((root/'examples/stages-seepage-column.json').read_text(encoding='utf-8'));mechanics['stages']=[];mechanics['verification_reject_once_stages']=[]
    flow=json.loads((root/'examples/richards-rain-pause.json').read_text(encoding='utf-8'))
    for m in flow['materials']:
        m['retention']['saturated_water_content']=.6 if mcc else .7/1.7
        m['retention']['saturated_specific_storage_per_m']=0.
        if mcc:m['saturated_conductivity_m_s']=[[1e-4,0.,0.],[0.,1e-4,0.],[0.,0.,1e-4]]
    flow['initial_head']={'head_m':1.5 if mcc else 0.}
    for i,p in enumerate(flow['analysis']['periods']):
        p['end_time_s']=(i+1)*1000.
        if mcc:
            p['name']=('drain','rain','drain_again')[i]
            p['boundaries'][1]['condition']['rainfall_vertical_m_s']=1e-5 if i==1 else 0.
    flow['settings']=dict(initial_time_step_s=100.,maximum_time_step_s=100.,minimum_time_step_s=.0001,store_integration_points=False)
    case=dict(schema='soil-fem-one-way/1',flow=flow,mechanics=mechanics,bishop=dict(kind='effective_saturation'),water_weight='water_equivalent_storage',solid_weight_reference='dry_bulk_unit_weight',porosity_policy='fixed_hydraulic_material',
        initial_porosity_tolerance=.02 if mcc else 1e-8,maximum_initial_engineering_strain=.05,initial_exterior_water=[],store_transfer_integration_points=True,
        exterior_periods=[dict(name=p['name'],reservoirs=[]) for p in flow['analysis']['periods']])
    return ground,case


if __name__=='__main__':
    for mcc,label in [(False,'mc'),(True,'mcc')]:
        ground,case=build(mcc)
        for p,data in [(root/f'examples/ground-one-way-{label}.json',ground),(root/f'examples/one-way-{label}.json',case)]:
            with p.open('x',encoding='utf-8') as f:json.dump(data,f,indent=2)
