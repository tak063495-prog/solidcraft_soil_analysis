"""New-only wet JOINT examples with explicit source/structure policies."""
import copy
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import fem_reference as fem

def make():
    def read(n):return json.loads((ROOT/'examples'/n).read_text(encoding='utf-8'))
    def write(n,v):
        with (ROOT/'examples'/n).open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
    for kind in ('mc','mcc'):
        g=read(f'ground-joint-{kind}.json');g['initialization']['water_unit_weight_kn_m3']=10.;g['initialization']['water_level_m']=.75;write(f'ground-joint-water-{kind}.json',g)
    for kind,source in (('mesh','joint-column-mesh.json'),('tet10-mesh','joint-column-tet10-mesh.json')):
        m=read(source)
        for axis,names in enumerate((('xmin','xmax'),('ymin','ymax'),('bottom','top'))):
            values=[min(x[axis] for x in m['nodes']),max(x[axis] for x in m['nodes'])]
            for side,n in enumerate(names):m['face_sets'][n]=[list(f) for f,owners in fem.topology(m).items() if len(owners)==1 and all(abs(m['nodes'][i][axis]-values[side])<1e-12 for i in f)]
        write(f'joint-water-column-{kind}.json',m)
    c=read('stages-joint-normal-column.json');c['verification_reject_once_stages']=[];c['joint_system']['joints'][0]['water']=dict(kind='continuous');c['stages']=[dict(kind='hydrostatic_water',name='rise',target_level_m=1.),dict(kind='pressure',name='preload',target_pressure_kpa=50.),dict(kind='hydrostatic_water',name='drawdown',target_level_m=.5)];write('stages-joint-water-column.json',c)
    materials=[dict(id=i,hydraulic_conductivity_m_s=[[k,0.,0.],[0.,k,0.],[0.,0.,k]]) for i,k in ((1,1e-4),(2,2e-4))]
    flow=dict(schema='soil-fem-darcy/1',units='m-s',materials=materials,heads=[dict(face_set='bottom',head_m=2.),dict(face_set='top',head_m=3.)])
    c['stages']=[dict(kind='steady_seepage',name='flow',flow=flow,exterior_water_level_m=.75),dict(kind='pressure',name='preload',target_pressure_kpa=20.),dict(kind='hydrostatic_water',name='hydrostatic_reset',target_level_m=1.)];write('stages-joint-water-flow.json',c)
    structure=copy.deepcopy(c);structure['joint_system']['structure_bodies']=['UpperClay'];structure['joint_system']['joints'][0]['water']=dict(kind='soil',soil_side='a');structure['stages'][0]['flow']['materials']=materials[:1];structure['stages'][0]['flow']['heads']=flow['heads'][:1];write('stages-joint-water-soil.json',structure)
    prescribed=copy.deepcopy(c);prescribed['joint_system']['joints'][0]['water']=dict(kind='prescribed',level_m=None);write('stages-joint-water-prescribed.json',prescribed)
    cycle=read('stages-joint-column.json');cycle['verification_reject_once_stages']=['shear','open','reverse'];cycle['joint_system']['joints'][0]['water']=dict(kind='continuous');write('stages-joint-water-cycle.json',cycle)

def make_local():
    def read(n):return json.loads((ROOT/'examples'/n).read_text(encoding='utf-8'))
    def write(n,v):
        with (ROOT/'examples'/n).open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
    design=read('local-excavation-cylinder.json');tool=design['commands'][-1]['params'];tool['base']=[5000,5000,3000];tool['height']=6000;write('joint-water-local-cylinder.json',design)
    g=read('ground-joint-local-mc.json');g['excavation_bodies']=['Pit'];g['initialization']['water_unit_weight_kn_m3']=10.;g['initialization']['water_level_m']=5.;write('ground-joint-water-local-mc.json',g)
    c=read('stages-joint-normal-column.json');c['monitors']={};c['monitor_sets']={'floor':'Pit:floor','surface':'Pit:surrounding_surface'};c['verification_reject_once_stages']=['wet_cut'];c['joint_system']['joints'][0]['water']=dict(kind='continuous');c['stages']=[dict(kind='hydrostatic_water',name='rise',target_level_m=6.),dict(kind='pressure',name='preload',target_pressure_kpa=30.),dict(kind='hydrostatic_water',name='drawdown',target_level_m=5.),dict(kind='excavate',name='wet_cut',bodies=['Pit'],remaining_surface_pressure_kpa=30.)];write('stages-joint-water-local.json',c)

if __name__=='__main__':
    if '--local-only' not in sys.argv:make()
    make_local()
