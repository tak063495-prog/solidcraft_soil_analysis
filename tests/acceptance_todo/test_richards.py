"""A-04 saturated equilibrium limit, not acceptance of general wet Richards JOINT.

Sources: interim/2026-11-16-issue-triage/todo.ja.md A-04,
interim/2026-11-13-joint-ext-combos/README.ja.md (a); src/richards.rs
Settings::default (head 1e-8 m, absolute flow 1e-16 m3/s, relative 1e-8).
The separate manufactured transient reference is in test_oracles.py.
"""
import copy
import json
import numpy as np
import pytest
from checks import CapabilityUnavailable, read, invoke, require_run
from verify_joint_ext_combos import wet_fixture
from joint_water_reference import Water


def richards_fixture(head, joint):
    mesh,_,wet=wet_fixture()
    if not joint:
        mesh=read('seepage-column-tet10-mesh.json')
    ground=read('ground-one-way-mc.json')
    case=read('one-way-mc.json')
    case['flow']['initial_head']={'head_m':head}
    case['flow']['analysis']['periods']=[dict(name='hold',end_time_s=1.,boundaries=[
        dict(face_set=face,condition=dict(kind='head',value=dict(head_m=head))) for face in ('bottom','top')])]
    case['flow']['settings']['initial_time_step_s']=.5
    case['flow']['settings']['maximum_time_step_s']=.5
    case['exterior_periods']=[dict(name='hold',reservoirs=[])]
    case['mechanics']['pressure_on_original_surface_only']=True
    case['mechanics']['stages']=[]
    assert ground['initialization']['water_level_m'] is None
    if joint:
        case['mechanics']['settings']['kinematics']='explicit'
        case['mechanics']['boundary']=copy.deepcopy(wet['boundary'])
        case['mechanics']['joint_system']=copy.deepcopy(wet['joint_system'])
        case['mechanics']['joint_system']['joints'][0].pop('birth_policy')
    return mesh,ground,case


def hydrostatic_richards_gates(mesh,ground,case,run,joint):
    assert run['complete'] is True
    head=case['flow']['initial_head']['head_m']
    frames=run['frames'];assert frames
    assert frames[-1]['state']['hydraulic']['time_s']==1.
    balances=[]
    for f in frames:
        state=f['state'];fluid=state['fluid']
        for key in ('hydraulic_points','mechanical_points'):
            points=fluid[key];assert points
            for p in points:
                expected=head-p['position_m'][2]
                assert abs(p['pressure_head_m']-expected)<=1e-8
                assert abs(p['physical_pore_pressure_kpa']-ground['initialization']['water_unit_weight_kn_m3']*expected)<=1e-8
                assert p['bishop_chi']==1. and p['effective_saturation']==1.
        balance=f.get('hydraulic_balance')
        if balance is not None:
            balances.append(balance)
            assert balance['maximum_scaled_residual']<=1.
            # Exact stationary manufactured limit: no storage or boundary flux.
            assert abs(balance['storage_rate_m3_s'])<=1e-16
            assert abs(balance['boundary_outflow_m3_s'])<=1e-16
            assert abs(balance['absolute_global_mass_residual_m3_s'])<=1e-16
        audit=f.get('mechanical_audit')
        if audit:
            assert audit['relative_effective_rhs_residual']<=1e-8
            assert audit['maximum_constraint_error_m']<=1e-10
            assert all(l['true_relative_residual']<=1e-10 for l in audit['linear_solves'])
        if joint:
            # Input-only independent face integration, never soil-only flow
            # success or an archived passed flag in place of joint mapping.
            water=Water(mesh,ground,case['mechanics'])
            mechanical=state['mechanical']
            expected=water.hydrostatic(mechanical['active_elements'],head)
            native=mechanical['joint_water']
            np.testing.assert_allclose(native['physical_pressure_kpa'],expected['pressure'],atol=1e-8,rtol=0.)
            np.testing.assert_allclose(native['nodal_load_kn'],expected['load'].ravel(),atol=1e-10,rtol=0.)
    assert balances, 'initialization alone is not a completed physical-time solve'


def test_A04_soil_only_saturated_limit_positive_control(solver,tmp_path):
    """Exercise the oracle and output schema without the missing JOINT coupling."""
    mesh,ground,case=richards_fixture(2.,False)
    result,output,_=invoke(solver,tmp_path,mesh,ground,case,'solve-one-way')
    assert result.returncode==0,result.stderr
    hydrostatic_richards_gates(mesh,ground,case,json.loads(output.read_text(encoding='utf-8')),False)


@pytest.mark.parametrize('head',[2.,3.])
@pytest.mark.xfail(strict=True,raises=CapabilityUnavailable,
    reason='A-04: ONE_WAY wet JOINT hydraulic mapping/coupled policy is unsupported; TODO A-04 / README(a), not soil-only Richards acceptance')
def test_A04_wet_joint_saturated_pressure_mapping_and_balance(solver,tmp_path,head):
    """todo A-04 / (a): independent head/pressure map and hydraulic balance.

    Constant total head 2/3 m, end time 1 s, head gate 1e-8 m, original
    absolute-flow gate 1e-16 m3/s and all-node pressure load gate 1e-10.
    """
    mesh,ground,case=richards_fixture(head,True)
    result,output,_=invoke(solver,tmp_path,mesh,ground,case,'solve-one-way')
    run=require_run(result,output,r'ONE_WAY requires dry-reference weights, explicit fluid/exterior loads and fixed active soil; construction/structures need their coupled policies')
    hydrostatic_richards_gates(mesh,ground,case,run,True)
