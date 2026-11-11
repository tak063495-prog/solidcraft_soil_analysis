"""A-01: two disconnected copies give an input-only 3D manufactured reference.

Sources: interim/2026-11-16-issue-triage/todo.ja.md A-01 and
interim/2026-11-13-joint-ext-combos/README.ja.md, independent elastic relative
1e-7, all-node loads 1e-10, actual restart exact except elapsed_seconds.
This probes multiple declared-tip JOINTs; it is not interacting-joint or wet
multiple-joint acceptance. The old single-JOINT verifier stays unchanged.
"""
import copy
import numpy as np
import pytest
import fem_reference as fem
from checks import CapabilityUnavailable,invoke,require_run,solver_gates,restart_exact
from verify_joint_ext_combos import tip_fixture
from joint_ext_reference import tip_elastic,elastic_operators


def two_copies(mesh,ground,case):
    mesh,ground,case=copy.deepcopy((mesh,ground,case))
    assert case['settings']['kinematics']=='explicit'
    assert all(not np.asarray(c.get('gradient',[0.,0.,0.])).any() for c in case['boundary']['constraints'])
    n=len(mesh['nodes']); e=len(fem.cells(mesh)); shift=np.array([4.,0.,0.])
    original=copy.deepcopy(mesh)
    mesh['nodes']+=(np.asarray(original['nodes'])+shift).tolist()
    for key in ('tetrahedra','tetrahedra10'):
        if key in original:mesh[key]+=(np.asarray(original[key],int)+n).tolist()
    mesh['material_ids']+=original['material_ids']
    region_shift=max(r['id'] for r in original['regions'])
    mesh['region_ids'] += [r+region_shift for r in original['region_ids']]
    for region in original['regions']:
        second=copy.deepcopy(region);second['id']+=region_shift;second['body']+=':copy2'
        second['bounds_m']=(np.asarray(second['bounds_m'])+shift).tolist()
        mesh['regions'].append(second)
    sides={case['joint_system']['joints'][0][s] for s in ('side_a','side_b')}
    for name,values in original['node_sets'].items():
        other=[i+n for i in values]
        if name in sides:mesh['node_sets'][name+':copy2']=other
        else:mesh['node_sets'][name]+=other
    for name,faces in original['face_sets'].items():
        other=(np.asarray(faces,int)+n).tolist()
        if name in sides:mesh['face_sets'][name+':copy2']=other
        else:mesh['face_sets'][name]+=other
    for body in copy.deepcopy(ground['bodies']):
        body['body']+=':copy2';ground['bodies'].append(body)
    states=ground['initialization']['geostatic']['seed']['states']
    for body,value in list(states.items()):states[body+':copy2']=copy.deepcopy(value)
    second=copy.deepcopy(case['joint_system']['joints'][0])
    second['id']='Seam2';second['side_a']+=':copy2';second['side_b']+=':copy2'
    tip=second['crack_tip']
    tip['edges']=(np.asarray(tip['edges'],int)+n).tolist()
    for face in tip['faces']:
        for side in ('a','b'):
            face['face_'+side]=[i+n for i in face['face_'+side]]
            face['owner_'+side]+=e
    case['joint_system']['joints'].append(second)
    for stage in case['stages']:
        if 'joints' in stage:stage['joints']=['Seam','Seam2']
    return mesh,ground,case


@pytest.mark.parametrize('policy',['transfer_force','stress_free'])
@pytest.mark.xfail(strict=True,raises=CapabilityUnavailable,
    reason='A-01: multiple declared-tip JOINT BIRTHs are restricted to one dry interface; TODO A-01 / README(a)')
def test_A01_multiple_tip_birth_input_only_full_fields_and_restart(solver,tmp_path,policy):
    """Two translated copies must reproduce the same independent elastic fields.

    (a): all 162 nodes, all volume GP stresses and live contact GP tractions
    relative 1e-7; each nodal physical load absolute 1e-10; real CLI restart.
    The reference is computed from the original input before reading the Run.
    """
    one,soil,single=tip_fixture(policy)
    expected=tip_elastic(one,soil,single)
    operators=elastic_operators(one,soil)
    bodies={r['id']:r['body'] for r in one['regions']}
    seed=soil['initialization']['geostatic']['seed']['states']
    mesh,ground,case=two_copies(one,soil,single)
    result,output,paths=invoke(solver,tmp_path,mesh,ground,case)
    run=require_run(result,output,r'declared-tip JOINT BIRTH is limited to one dry interface')
    solver_gates(run)
    assert len(run['frames'])==len(expected)
    for frame,reference in zip(run['frames'],expected,strict=True):
        u=np.tile(reference['u'].reshape(-1,3),(2,1))
        assert np.linalg.norm(np.asarray(frame['state']['displacement_m'])-u)/max(np.linalg.norm(u),1e-30)<=1e-7
        stress=[]
        for i,(_,d,ids,_) in enumerate(operators):
            q=fem.volume(one,i)
            s0=seed[bodies[one['region_ids'][i]]]['effective_stress_compression_positive_kpa']
            stress.extend(np.asarray(s0)-np.einsum('ab,qbi,i->qa',d,q['b'],reference['u'][ids]))
        stress=np.tile(stress,(2,1))
        actual=np.asarray([p['effective_stress_compression_positive_kpa'] for p in frame['state']['material_states']])
        assert np.linalg.norm(actual-stress)/max(np.linalg.norm(stress),1e-30)<=1e-7
        np.testing.assert_allclose(frame['state']['applied_nodal_load_kn'],np.tile(reference['load'],2),atol=1e-10,rtol=0.)
        if reference['born']:
            traction=np.tile(reference['traction'],(2,1))
            actual=np.asarray([p['traction_kpa'] for p in frame['state']['joints']['points']])
            assert np.linalg.norm(actual-traction)/max(np.linalg.norm(traction),1.)<=1e-7
    restart_exact(solver,tmp_path,paths,run)
