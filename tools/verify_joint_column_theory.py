"""A contact-supported upper layer: constrained-modulus + finite-Kn theory.

Upper-layer vertical translation is free; only JOINT supports that mode.
Both layers have lateral rollers, supplied either directly or by explicit
node-set constraints. This does not test a block free in all six rigid modes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def verify(run,output):
    mp=run/'ground/mesh.json';gp=run/'ground/ground.json';rp=run/'stage-results.json'
    mesh=json.loads(mp.read_text());ground=json.loads(gp.read_text());native=json.loads(rp.read_text());case=native['case']
    assert native['complete']
    if case['settings']['kinematics']=='explicit':
        bc=case['boundary'];assert not bc.get('nodal_loads') and not bc.get('surface_loads')
        expected={('xmin',(1,0,0)),('xmax',(1,0,0)),('ymin',(0,1,0)),('ymax',(0,1,0)),('bottom',(0,0,1))}
        assert {(c['node_set'],tuple(c['direction'])) for c in bc['constraints']}==expected
        assert len(bc['constraints'])==5 and all(c['value_m']==0 and not any(c.get('gradient',[0,0,0])) for c in bc['constraints'])
    else:
        assert case['settings']['kinematics']=='side_rollers' and not case.get('boundary')
    xyz=np.asarray(mesh['nodes']);assert xyz[:,2].min()==0 and xyz[:,2].max()==1
    materials=ground['materials'];assert all(m['mohr_coulomb']==materials[0]['mohr_coulomb'] for m in materials)
    assert all(m['unit_weight_kn_m3']==20 for m in materials)
    e=materials[0]['mohr_coulomb']['young_kpa'];nu=materials[0]['mohr_coulomb']['poisson'];modulus=e*(1-nu)/((1+nu)*(1-2*nu))
    kn=case['joint_system']['joints'][0]['material']['normal_stiffness_kpa_m'];regions={r['id']:r['body'] for r in mesh['regions']}
    upper=set(n for ei,t in enumerate(ref.cells(mesh)) if regions[mesh['region_ids'][ei]]=='UpperClay' for n in t)
    side=np.asarray([n in upper for n in range(len(xyz))]);z=xyz[:,2];rows=[]
    for frame in native['frames']:
        state=frame['state'];active=np.asarray(state['active_elements']);used=np.unique(ref.cells(mesh)[active]);p=state['surface_pressure_kpa']
        cut=frame['name']=='excavate_upper';uz=((90*z+10*z*z)/modulus if cut else -(p*z+20*(z-z*z/2))/modulus-side*(110+p)/kn)
        u=np.zeros_like(xyz);u[:,2]=uz;actual=np.asarray(state['displacement_m']);error=np.linalg.norm(actual[used]-u[used])/np.linalg.norm(u[used])
        assert error<1e-7,(frame['name'],error)
        stress_error=0.;void_error=0.
        for point,s in zip(native['integration_points'],state['material_states']):
            if not active[point['element_index']]:continue
            height=point['position_m'][2];sz=20*(.5-height) if cut else 100+p+20*(1-height)
            expected=[100+nu/(1-nu)*(sz-100)]*2+[sz,0.,0.,0.]
            stress_error=max(stress_error,float(np.max(np.abs(np.asarray(s['effective_stress_compression_positive_kpa'])-expected))))
            void_error=max(void_error,abs(s['void_ratio']-(1.7*np.exp(-(sz-100)/modulus)-1)))
            assert s['equivalent_plastic_strain']==0
        assert stress_error<1e-5 and void_error<1e-8,(stress_error,void_error)
        if cut:
            assert state['joints']==native['frames'][-2]['state']['joints']
            assert all(not p['active'] for p in frame['joint_results']['points'])
            holding=np.asarray(native['excavation_audits'][0]['initial_holding_load_kn']).reshape(-1,3).sum(axis=0)
            assert np.allclose(holding,[0,0,-110],atol=1e-6)
        else:
            for point in frame['joint_results']['points']:
                assert point['active'] and abs(point['normal_pressure_kpa']-(110+p))<1e-6
                assert abs(point['gap_m']+(110+p)/kn)<1e-10
        rows.append({'frame':frame['name'],'relative_all_active_node_3d_displacement_error':float(error),
            'maximum_stress_error_kpa':stress_error,'maximum_void_ratio_error':void_error})
    result={'schema':'soil-fem-joint-normal-column-theory/1','reference':'analytic constrained-modulus soil strain, finite normal interface stiffness, exact depth-dependent self-weight stress; free upper vertical mode supported by JOINT',
        'frames':rows,'solver_sha256':native['provenance']['solver_binary_sha256'],
        'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (mp,gp,rp,Path(__file__))}}
    with output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);a=p.parse_args();verify(a.run.resolve(),a.output.resolve())
