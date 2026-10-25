"""New-only fixtures for initially inactive fill and explicit MATSET water policy."""
import copy
import json
from pathlib import Path

def make(root):
    ground=json.loads((root/'examples/ground-one-way-mc.json').read_text(encoding='utf-8'))
    extra=copy.deepcopy(ground['materials'][1]);extra['id']=3;extra['name']='Improved fill';extra['mohr_coulomb']['young_kpa']=20000.;ground['materials'].append(extra);ground['construction_material_ids']=[3]
    base=json.loads((root/'examples/one-way-mc.json').read_text(encoding='utf-8'))
    material=copy.deepcopy(base['flow']['materials'][1]);material['id']=3;material['retention']['saturated_water_content']=.5;material['saturated_conductivity_m_s']=[[5e-6,0,0],[0,5e-6,0],[0,0,5e-6]];base['flow']['materials'].append(material)
    base['initial_porosity_tolerance']=.15
    periods=[]
    for i,name in enumerate(('hold','born','changed')):
        boundaries=[dict(face_set='bottom',condition=dict(kind='head',value=dict(head_m=0.))),dict(face_set='floor' if i==0 else 'top',condition=dict(kind='atmospheric_no_ponding',rainfall_vertical_m_s=0.))]
        boundaries += [dict(face_set=n,condition=dict(kind='outward_flux',flux_m_s=0.)) for n in ('xmin','xmax','ymin','ymax')]
        periods.append(dict(name=name,end_time_s=(i+1)*100.,boundaries=boundaries))
    base['flow']['analysis']=dict(kind='transient',periods=periods);base['exterior_periods']=[dict(name=p['name'],reservoirs=[]) for p in periods]
    fresh=dict(effective_stress_compression_positive_kpa=[0.]*6,plastic_strain_compression_positive=[0.]*6,void_ratio=.7,equivalent_plastic_strain=0.,plastic_work_kpa=0.,last_return='elastic')
    common=dict(water_boundary_policy='explicit_period_boundaries',surface_traction_policy='explicit_period_reservoirs')
    case=dict(schema='soil-fem-one-way-construction/1',one_way=base,excavations=[],initially_inactive_bodies=['UpperClay'],transitions=[dict(common,name='place',before_period='born',action=dict(kind='birth',bodies=['UpperClay'],states=dict(UpperClay=fresh),placement='cad_new_nodes',initial_head=dict(head_m=0.))),dict(common,name='improve',before_period='changed',action=dict(kind='material_set',bodies=['UpperClay'],material_id=3,history=dict(mode='preserve'),water=dict(kind='preserve_water',reference_head=dict(head_m=0.))))])
    for name,value in (('ground-one-way-lifecycle',ground),('one-way-lifecycle',case)):
        with (root/f'examples/{name}.json').open('x',encoding='utf-8') as stream:json.dump(value,stream,indent=2)

if __name__=='__main__':make(Path(__file__).resolve().parents[1])
