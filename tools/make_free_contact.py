"""Input-only free upper block on a horizontal or inclined dry interface."""
import argparse
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def make(angle=0.,friction=.3,mcc=False,force=(0.,0.,0.),tangential_stiffness=None,pc=None,preload_kpa=5.):
    root=Path(__file__).resolve().parents[1]
    mesh=json.loads((root/'examples/joint-column-tet10-mesh.json').read_text())
    ground=json.loads((root/'examples'/('ground-joint-mcc.json' if mcc else 'ground-joint-mc.json')).read_text())
    case=json.loads((root/'examples/stages-joint-explicit-normal-column.json').read_text())
    xyz=np.asarray(mesh['nodes']);top=[list(f) for f,o in ref.topology(mesh).items() if len(o)==1 and np.all(abs(xyz[list(f),2]-1)<1e-10)]
    mesh['face_sets']['top']=top
    a=np.deg2rad(angle);rotation=np.array([[np.cos(a),0,np.sin(a)],[0,1,0],[-np.sin(a),0,np.cos(a)]])
    xyz=xyz@rotation.T;mesh['nodes']=xyz.tolist()
    for region in mesh['regions']:
        nodes=np.unique(ref.cells(mesh)[np.asarray(mesh['region_ids'])==region['id']])
        region['bounds_m']=[xyz[nodes].min(axis=0).tolist(),xyz[nodes].max(axis=0).tolist()]
    case['boundary']={'constraints':[{'node_set':'bottom','direction':direction,'value_m':0.} for direction in np.eye(3).tolist()],
        'surface_loads':[{'kind':'traction','face_set':'top','traction_kpa':[0.,0.,-preload_kpa]}]}
    case['joint_system']['joints'][0]['material']['friction_coefficient']=friction
    if tangential_stiffness is not None:case['joint_system']['joints'][0]['material']['tangential_stiffness_kpa_m']=tangential_stiffness
    case['settings']['maximum_newton_iterations']=60
    case['stages']=[{'kind':'boundary','name':'load','target_factor':1.,'initial_increment':.5},
                    {'kind':'boundary','name':'unload','target_factor':0.,'initial_increment':.5},
                    {'kind':'excavate','name':'excavate_upper','bodies':['UpperClay'],'remaining_surface_pressure_kpa':0.,'initial_increment':.5}]
    case['verification_reject_once_stages']=['load','excavate_upper']
    geo=ground['initialization']['geostatic'];geo['verification_reject_once']=True
    if not mcc:
        geo['base_loads']={}
        for state in geo['seed']['states'].values():state['effective_stress_compression_positive_kpa']=[0.]*6
    if any(force):
        area=sum(np.linalg.norm(np.cross(xyz[f[1]]-xyz[f[0]],xyz[f[2]]-xyz[f[0]]))/2 for f in top)
        geo['base_loads'].setdefault('surface_loads',[]).append({'kind':'traction','face_set':'top','traction_kpa':(np.asarray(force)/area).tolist()})
    if pc is not None:
        if not mcc or not np.isfinite(pc) or pc<100:raise ValueError('MCC seed requires finite pc >= 100 kPa')
        materials={m['id']:m['mcc'] for m in ground['materials']};bodies={b['body']:materials[b['material_id']] for b in ground['bodies']}
        for name,state in geo['seed']['states'].items():
            p=bodies[name];v=p['specific_volume_ncl_at_reference']-p['lambda']*np.log(pc/p['reference_pressure_kpa'])+p['kappa']*np.log(pc/100)
            if v<=1:raise ValueError('MCC seed void ratio must be positive')
            state['preconsolidation_pressure_kpa']=pc;state['void_ratio']=float(v-1)
    return mesh,ground,case


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--angle',type=float,default=0.);p.add_argument('--friction',type=float,default=.3);p.add_argument('--mcc',action='store_true');p.add_argument('--force',type=float,nargs=3,default=(0.,0.,0.));p.add_argument('--tangential-stiffness',type=float);p.add_argument('--pc',type=float);p.add_argument('--preload-kpa',type=float,default=5.)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    for name,value in zip(('mesh','ground','stages'),make(a.angle,a.friction,a.mcc,a.force,a.tangential_stiffness,a.pc,a.preload_kpa)):
        (a.output/(name+'.json')).write_text(json.dumps(value,indent=2),encoding='utf-8')
