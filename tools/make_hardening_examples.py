"""Explicit example inputs for IDH; no material parameters are inferred from CAD."""
import copy
import json
from pathlib import Path


def initial(stress,e,pc):
    return {'effective_stress_compression_positive_kpa':stress,'plastic_strain_compression_positive':[0.]*6,
        'void_ratio':e,'shear_plastic_strain':0.,'cap_plastic_volumetric_strain':0.,
        'initial_cap_pressure_kpa':pc,'cap_pressure_kpa':pc,'plastic_work_kpa':0.,'last_return':'elastic'}


def make():
    root=Path(__file__).resolve().parents[1];examples=root/'examples'
    read=lambda n:json.loads((examples/n).read_text())
    p=read('hardening-parameters.json');ground=read('ground-mc-column.json');case=read('stages-mc-column.json')
    for i,material in enumerate(ground['materials']):
        material.pop('mohr_coulomb');material['name']=f'IDH example layer {i+1}'
        material['unit_weight_kn_m3']=20.;material['invariant_double_hardening']=copy.deepcopy(p)
        if i==1:
            material['invariant_double_hardening']['unloading_young_ref_kpa']=50000.
            material['invariant_double_hardening']['cap_hardening_modulus_ref_kpa']=8000.
    for layer in ground['layers']:
        layer['k0_x']=layer['k0_y']=1.
        layer['consolidation']={'definition':'hardening_cap_over_mean_effective','ratio':1.05}
    case.pop('boundary');case['settings']['kinematics']='oedometer';case['settings']['maximum_newton_iterations']=40
    case['settings']['material_maximum_substep']=.0002
    case['monitors']={'ground_surface':1.,'excavation_floor':.5}
    case['stages']=[{'kind':'pressure','name':'preload','target_pressure_kpa':200.},
        {'kind':'pressure','name':'unload','target_pressure_kpa':100.},
        {'kind':'pressure','name':'reload','target_pressure_kpa':150.},
        {'kind':'excavate','name':'excavate_upper','bodies':['UpperClay'],'remaining_surface_pressure_kpa':100.,'initial_increment':1.}]
    case['verification_reject_once_stages']=['preload','excavate_upper']
    mesh=read('mc-column-mesh.json');mesh['schema']='soil-fem-mesh/3';mesh['quadrature_order']=3
    edges={};mesh['tetrahedra10']=[]
    for cell in mesh['tetrahedra']:
        mids=[]
        for a,b in ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3)):
            edge=tuple(sorted((cell[a],cell[b])))
            if edge not in edges:
                edges[edge]=len(mesh['nodes']);mesh['nodes'].append([(x+y)/2 for x,y in zip(mesh['nodes'][edge[0]],mesh['nodes'][edge[1]])])
            mids.append(edges[edge])
        mesh['tetrahedra10'].append(cell+mids)
    mesh['node_sets']={name:[i for i,p in enumerate(mesh['nodes']) if (name=='all' or (name=='top' and p[2]==1.) or (name=='bottom' and p[2]==0.))] for name in ('all','top','bottom')}
    path={'schema':'soil-fem-hardening-path/1','parameters':p,'initial_state':initial([100.,100.,100.,0.,0.,0.],.7,120.),
        'maximum_substep':.0001,'increments':[[.001,.001,.001,0.,0.,0.],[.004,-.002,-.002,.001,.0005,0.],[-.0005,-.0005,-.0005,0.,0.,0.]]}
    local=read('ground-local.json');local_case=read('stages-local-sparse.json')
    for material in local['materials']:
        material.pop('mcc');material['name']='Illustrative cohesive IDH layer'
        material['invariant_double_hardening']={**p,'shear_pressure_shift_kpa':30.}
    for layer in local['layers']:
        layer['consolidation']={'definition':'hardening_cap_over_mean_effective','ratio':1.5}
        layer['void_ratio']={'mode':'explicit','value':.7}
    local_case['settings']['maximum_newton_iterations']=40
    local_case['settings']['material_maximum_substep']=.0002
    local_case['stages']=[{'kind':'pressure','name':'preload','target_pressure_kpa':90.,'initial_increment':.5},
        {'kind':'pressure','name':'unload','target_pressure_kpa':30.,'initial_increment':.5},local_case['stages'][0]]
    local_case['verification_reject_once_stages']=['preload','excavate_pit']
    return {'ground-hardening-column.json':ground,'stages-hardening-column.json':case,
        'hardening-column-tet10-mesh.json':mesh,'hardening-path.json':path,
        'ground-hardening-local.json':local,'stages-hardening-local.json':local_case}


if __name__=='__main__':
    destination=Path(__file__).resolve().parents[1]/'examples'
    for name,value in make().items():
        path=destination/name
        if path.exists():assert json.loads(path.read_text())==value,'refuse to overwrite a changed example'
        else:
            with path.open('x',encoding='utf-8') as stream:json.dump(value,stream,indent=2)
