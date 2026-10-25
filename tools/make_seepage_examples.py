"""Deterministic confined columns with an exposed post-excavation floor."""
import copy
import json
from pathlib import Path
import fem_reference as ref


def tagged(mesh):
    xyz=mesh['nodes'];faces=ref.topology(mesh)
    sets={}
    for axis,lo,hi in ((0,'xmin','xmax'),(1,'ymin','ymax'),(2,'bottom','top')):
        for level,name in ((0.,lo),(1.,hi)):
            sets[name]=[list(f) for f,o in faces.items() if len(o)==1 and all(abs(xyz[n][axis]-level)<1e-12 for n in f)]
    sets['floor']=[list(f) for f,o in faces.items() if len(o)==2 and all(abs(xyz[n][2]-.5)<1e-12 for n in f)]
    mesh['face_sets']=sets
    return mesh


def flow(bottom=1.4,top=1.,surface='top'):
    return {'schema':'soil-fem-darcy/1','units':'m-s','materials':[
        {'id':i,'hydraulic_conductivity_m_s':[[1e-6,0.,0.],[0.,1e-6,0.],[0.,0.,1e-6]]} for i in (1,2)],
        'heads':[{'face_set':'bottom','head_m':bottom},{'face_set':surface,'head_m':top}]}


def stages():
    root=Path(__file__).resolve().parents[1]
    case=json.loads((root/'examples/stages-water-column.json').read_text())
    case['stages']=[
        {'kind':'steady_seepage','name':'upward_flow','flow':flow(),'exterior_water_level_m':1.,'initial_increment':.5},
        {'kind':'steady_seepage','name':'downward_flow','flow':flow(.6),'exterior_water_level_m':1.,'initial_increment':1.},
        {'kind':'hydrostatic_water','name':'recover_hydrostatic','target_level_m':1.,'initial_increment':.5},
        {'kind':'excavate','name':'excavate_upper','bodies':['UpperClay'],'remaining_surface_pressure_kpa':100.,'initial_increment':.5},
        {'kind':'steady_seepage','name':'dewater_remaining','flow':flow(.1,.5,'floor'),'exterior_water_level_m':.5,'initial_increment':.5}]
    case['verification_reject_once_stages']=['upward_flow','downward_flow','excavate_upper','dewater_remaining']
    return case


def local_stages():
    root=Path(__file__).resolve().parents[1];case=json.loads((root/'examples/stages-local-sparse.json').read_text())
    flow_case=flow(8.5,8.)
    for material in flow_case['materials']:
        material['hydraulic_conductivity_m_s']=([[1e-6,.2e-6,0.],[.2e-6,2e-6,0.],[0.,0.,1e-6]] if material['id']==1 else [[2e-6,.1e-6,0.],[.1e-6,1e-6,0.],[0.,0.,3e-6]])
    after=copy.deepcopy(flow_case);after['heads']=[{'face_set':name,'head_m':8.5 if name=='bottom' else 8.} for name in ('bottom','Pit:surrounding_surface','Pit:floor','Pit:wall')]
    case['stages']=[
        {'kind':'pressure','name':'preload','target_pressure_kpa':90.,'initial_increment':1.},
        {'kind':'pressure','name':'unload','target_pressure_kpa':30.,'initial_increment':.5},
        {'kind':'steady_seepage','name':'recharge','flow':flow_case,'exterior_water_level_m':8.,'initial_increment':.5},
        {'kind':'excavate','name':'excavate_pit','bodies':['Pit'],'remaining_surface_pressure_kpa':30.,'initial_increment':.5},
        {'kind':'steady_seepage','name':'reflow_after_excavation','flow':after,'exterior_water_level_m':8.,'initial_increment':.5}]
    case['verification_reject_once_stages']=['preload','recharge','excavate_pit','reflow_after_excavation']
    return case


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    for name,source in (('tet4','mc-column-mesh.json'),('tet10','hardening-column-tet10-mesh.json')):
        mesh=tagged(json.loads((root/'examples'/source).read_text()))
        (root/f'examples/seepage-column-{name}-mesh.json').write_text(json.dumps(mesh,indent=2),encoding='utf-8')
    (root/'examples/stages-seepage-column.json').write_text(json.dumps(stages(),indent=2),encoding='utf-8')
    (root/'examples/flow-seepage-column.json').write_text(json.dumps(flow(),indent=2),encoding='utf-8')
    ground=json.loads((root/'examples/ground-water-column-mcc.json').read_text())
    for layer in ground['layers']:layer['consolidation']['ratio']=1.001
    (root/'examples/ground-seepage-column-mcc.json').write_text(json.dumps(ground,indent=2),encoding='utf-8')
    for kind in ('mc','mcc'):
        ground=json.loads((root/('examples/ground-water-local-mc.json' if kind=='mc' else 'examples/ground-water-local.json')).read_text())
        ground['initialization']['water_level_m']=8.
        for layer in ground['layers']:layer['k0_x']=layer['k0_y']=1.
        (root/f'examples/ground-seepage-local-{kind}.json').write_text(json.dumps(ground,indent=2),encoding='utf-8')
    (root/'examples/stages-seepage-local.json').write_text(json.dumps(local_stages(),indent=2),encoding='utf-8')
