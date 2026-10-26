"""Small two-layer K0/MCC/member construction specimen (no CAD acceptance)."""
import argparse
import copy
import json
from pathlib import Path


def inputs(root, projected=True, increment=.125):
    mesh=json.loads((root/'examples/hardening-column-tet10-mesh.json').read_text())
    ground=json.loads((root/'examples/ground-mixed-mcc.json').read_text())
    mesh['node_sets']['member_base']=[4]
    mesh['node_sets']['member_end']=[7]
    case=json.loads((root/'examples/stages-geostatic-column.json').read_text())
    case['schema']='soil-fem-projected-mcc-stages/1' if projected else 'soil-fem-mixed-mcc-stages/1'
    case['settings'].update(kinematics='side_rollers',initial_increment=increment,
        material_maximum_substep=.0001,minimum_increment=1/4096,maximum_newton_iterations=45)
    case['boundary']={'constraints':[]}
    case['pressure_on_original_surface_only']=False
    case['monitor_sets']={'surface':'top'}
    case['monitors']={'mid_height':.5,'top':1.}
    case['beam_system']={
        'schema':'soil-fem-structures/1',
        'nodes':[
            {'id':'a','position_m':mesh['nodes'][4],'soil_node_set':'member_base'},
            {'id':'b','position_m':mesh['nodes'][7],'soil_node_set':'member_end'},
            {'id':'c','position_m':[-1.,0.,.5],'soil_node_set':None},
            {'id':'d','position_m':[-1.,1.,.5],'soil_node_set':None}],
        'sections':[{'id':'strut','young_kpa':2e7,'poisson':.2,'area_m2':.005,
            'iy_m4':1e-5,'iz_m4':1e-5,'torsion_j_m4':1e-5,'unit_weight_kn_m3':24.,
            'formulation':{'kind':'euler_bernoulli'}}],
        'members':[{'id':'beam','nodes':['a','b'],'section_id':'strut','local_y':[0.,0.,1.],
            'initially_active':False,'boundary_line_load_kn_m':[0.,0.,-1.]}],
        'plate_sections':[{'id':'slab','young_kpa':3e7,'poisson':.2,'thickness_m':.1,
            'unit_weight_kn_m3':24.,'drilling_factor':1e-4}],
        'plates':[{'id':'plate','nodes':['a','c','d'],'section_id':'slab',
            'initially_active':False,'boundary_traction_kpa':[0.,0.,-2.]}],
        'rotation_length_scale_m':1.,'maximum_condensed_dofs':60,
        'rotation_constraints':[{'node_id':'a','direction':[int(i==j) for i in range(3)],'value':0.} for j in range(3)]}
    case['stages']=[
        {'kind':'pressure','name':'preload','target_pressure_kpa':140.},
        {'kind':'install_structures','name':'install','members':['beam','plate'],'placement':'cad_new_nodes'},
        {'kind':'boundary','name':'member_load','target_factor':1.},
        {'kind':'pressure','name':'hardening','target_pressure_kpa':180.},
        {'kind':'pressure','name':'unload','target_pressure_kpa':100.},
        {'kind':'excavate','name':'excavate_upper','bodies':['UpperClay'],'remaining_surface_pressure_kpa':100.},
        {'kind':'boundary','name':'member_unload','target_factor':0.},
        {'kind':'remove_structures','name':'remove','members':['beam','plate']}]
    case['verification_reject_once_stages']=['install','hardening','excavate_upper','remove']
    return copy.deepcopy((mesh,ground,case))


def prepare(root, output, projected=True, increment=.125):
    output.mkdir(parents=True,exist_ok=False)
    ground_dir=output/'ground';ground_dir.mkdir()
    for path,value in zip((ground_dir/'mesh.json',ground_dir/'ground.json',output/'stages.json'),
                          inputs(root,projected,increment),strict=True):
        with path.open('x',encoding='utf-8') as stream:json.dump(value,stream,indent=2,allow_nan=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path)
    p.add_argument('--gp-volume',action='store_true');p.add_argument('--increment',type=float,default=.125)
    a=p.parse_args();prepare(Path(__file__).resolve().parents[1],a.output.resolve(),not a.gp_volume,a.increment)
