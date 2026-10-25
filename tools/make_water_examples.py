"""Deterministic staged-water examples; retain existing different inputs."""
import copy
import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
def read(name):return json.loads((ROOT/'examples'/name).read_text())
def write(name,value):
    p=ROOT/'examples'/name
    if p.exists():assert json.loads(p.read_text())==value,p
    else:
        with p.open('x') as f:json.dump(value,f,indent=2)


def make():
    g=read('ground-mc-column.json');g['initialization'].update(water_level_m=1.,water_unit_weight_kn_m3=10.)
    for m in g['materials']:
        m['unit_weight_kn_m3']=20.;m['mohr_coulomb'].update(young_kpa=10000.,poisson=.25,cohesion_kpa=500.,dilation_deg=0.)
    for layer in g['layers']:layer.update(k0_x=1.,k0_y=1.)
    write('ground-water-column-mc.json',g)
    c=read('stages-hardening-column.json');c['monitor_sets']={};c['monitors']={'top':1.,'floor':.5};c['pressure_on_original_surface_only']=False
    c['settings']['initial_increment']=.5
    c['stages']=[{'kind':'hydrostatic_water','name':name,'target_level_m':h,'initial_increment':.5} for name,h in [('drawdown',.37),('flood',1.3),('drain',-.2),('recover',1.)]]
    c['stages'].append({'kind':'excavate','name':'excavate_upper','bodies':['UpperClay'],'remaining_surface_pressure_kpa':100.,'initial_increment':.5})
    c['verification_reject_once_stages']=['drawdown','excavate_upper']
    write('stages-water-column.json',c)
    mcc=read('ground.json');mcc['initialization'].update(surface_surcharge_kpa=100.,water_level_m=1.,water_unit_weight_kn_m3=10.)
    for m in mcc['materials']:m['unit_weight_kn_m3']=20.
    for layer in mcc['layers']:
        layer.update(k0_x=1.,k0_y=1.,consolidation={'definition':'mcc_pc_over_mean_effective','ratio':1.02})
    write('ground-water-column-mcc.json',mcc)
    local=read('ground-local.json');local['initialization']['water_unit_weight_kn_m3']=10.
    write('ground-water-local.json',local)
    stages=read('stages-local-sparse.json');stages['stages']=[
        {'kind':'hydrostatic_water','name':'raise_water','target_level_m':7.,'initial_increment':.5},
        stages['stages'][0],
        {'kind':'hydrostatic_water','name':'drawdown','target_level_m':6.7,'initial_increment':.5},
        {'kind':'hydrostatic_water','name':'reflood','target_level_m':7.3,'initial_increment':.5}]
    stages['verification_reject_once_stages']=['raise_water','excavate_pit','drawdown']
    write('stages-water-local.json',stages)
    plastic=copy.deepcopy(stages)
    plastic['stages']=[{'kind':'pressure','name':'preload','target_pressure_kpa':90.,'initial_increment':1.},
        {'kind':'pressure','name':'unload','target_pressure_kpa':30.,'initial_increment':.5}]+plastic['stages']
    plastic['verification_reject_once_stages'].insert(0,'preload')
    write('stages-water-local-mcc.json',plastic)
    elastic=read('ground-mc-local.json');elastic['initialization']['water_unit_weight_kn_m3']=10.
    for m in elastic['materials']:
        m['mohr_coulomb'].update(cohesion_kpa=500.,dilation_deg=0.,tension_cutoff_kpa=None)
    write('ground-water-local-mc.json',elastic)


if __name__=='__main__':make()
