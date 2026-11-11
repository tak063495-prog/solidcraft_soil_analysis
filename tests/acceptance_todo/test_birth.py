"""A-01/A-02/A-03, from interim/2026-11-16-issue-triage/todo.ja.md.

(a) = interim/2026-11-13-joint-ext-combos/README.ja.md,
sections "受理する範囲と残す拒否", "独立参照と保存再開".

Strict xfail only covers the named existing scope guard. Completion, residuals,
all-node water loads, work and actual restart checks follow successful solves.
These criteria are not a declaration that the entire TODO item is complete.
"""
import copy
import numpy as np
import pytest
from checks import (CapabilityUnavailable, read, invoke, require_run, refused,
                    solver_gates, contact_water_gates, restart_exact,
                    rehashed_birth_tampering_refused)
from verify_joint_ext_combos import tip_fixture, wet_fixture


@pytest.mark.parametrize('policy', ['transfer_force','stress_free'])
@pytest.mark.xfail(strict=True,raises=CapabilityUnavailable,
    reason='A-01: declared-tip wet JOINT BIRTH is still restricted to one dry interface; TODO A-01 / README(a)')
def test_A01_wet_tip_birth_nodal_water_work_and_real_restart(solver,tmp_path,policy):
    """(a): water/load reassembly 1e-10, head/pressure 1e-8, work 1e-10."""
    mesh,ground,case=tip_fixture(policy)
    ground['initialization']['water_level_m']=1.5
    case['joint_system']['joints'][0]['water']={'kind':'continuous'}
    result,output,paths=invoke(solver,tmp_path,mesh,ground,case)
    run=require_run(result,output,r'declared-tip JOINT BIRTH is limited to one dry interface')
    solver_gates(run)
    contact_water_gates(mesh,ground,case,run)
    restart_exact(solver,tmp_path,paths,run)
    rehashed_birth_tampering_refused(solver,tmp_path,paths,run)


@pytest.mark.parametrize('mode', ['crossing_body_birth','partly_active_contact'])
def test_A02_undefined_reactivation_remains_refused(solver,tmp_path,mode):
    """todo A-02 / (a): crossing and incomplete owners stay refused.

    Successful acceptance semantics are unspecified, so do not mark an unsafe
    activation as a desired pass or xfail. Assert the specific current refusal.
    """
    mesh,ground,case=tip_fixture('stress_free')
    case['initially_inactive_bodies']=['UpperClay']
    case['verification_reject_once_stages']=[]
    if mode=='crossing_body_birth':
        seed=copy.deepcopy(ground['initialization']['geostatic']['seed']['states']['UpperClay'])
        seed['effective_stress_compression_positive_kpa']=[0.]*6
        case['stages']=[dict(kind='birth',name='cross_tip',bodies=['UpperClay'],states={'UpperClay':seed},placement='cad_new_nodes')]
        pattern=r'body BIRTH across a declared crack tip is unverified'
    else:
        case['stages']=[dict(kind='joint_birth',name='incomplete_owners',joints=['Seam'])]
        pattern=r'JOINT BIRTH needs inactive contact and both active bodies'
    result,output,_=invoke(solver,tmp_path,mesh,ground,case)
    refused(result,output,pattern)


def test_A03_supported_hydrostatic_birth_audit_positive_control(solver,tmp_path):
    """(a) limited hydrostatic MC birth exercises the audit beyond xfail guards."""
    mesh,ground,case=wet_fixture()
    case['stages']=case['stages'][:3]
    case['verification_reject_once_stages']=[]
    result,output,paths=invoke(solver,tmp_path,mesh,ground,case)
    assert result.returncode==0,result.stderr
    import json
    run=json.loads(output.read_text(encoding='utf-8'))
    solver_gates(run)
    contact_water_gates(mesh,ground,case,run)
    restart_exact(solver,tmp_path,paths,run)
    rehashed_birth_tampering_refused(solver,tmp_path,paths,run)


@pytest.mark.parametrize('policy', ['continuous','soil','prescribed'])
@pytest.mark.xfail(strict=True,raises=CapabilityUnavailable,
    reason='A-03: steady Darcy plus soil BIRTH is unverified; TODO A-03 / README(a), no pressure/history auto-fill')
def test_A03_darcy_birth_hydrostatic_manufactured_limit(solver,tmp_path,policy):
    """Darcy H=1.5 => zero flux, p=gamma(1.5-z); original gates unchanged.

    This well-defined zero-gradient limit is deliberately separate from a
    nonzero-gradient birth initialization, whose new input contract is absent.
    """
    mesh,ground,case=wet_fixture()
    case['stages']=case['stages'][:3]
    case['verification_reject_once_stages']=[]
    definition=case['joint_system']['joints'][0]
    definition['water']={'kind':policy}
    if policy=='soil':
        # SOIL is a soil/impermeable-structure policy, never soil/soil.
        definition['water']['soil_side']='b'
        case['joint_system']['structure_bodies']=['LowerClay']
    if policy=='prescribed': definition['water']['level_m']=1.5
    flow=read('stages-joint-water-flow.json')['stages'][0]['flow']
    flow['heads']=[dict(face_set='bottom',head_m=1.5)]
    if policy=='soil': flow['heads']=[dict(face_set='top',head_m=1.5)]
    # Before birth the SOIL-side domain is empty. Solve Darcy after creating
    # it; the other policies have a pre-existing hydraulic soil domain.
    index=1 if policy=='soil' else 0
    case['stages'].insert(index,dict(kind='steady_seepage',name='darcy',flow=flow,exterior_water_level_m=1.5))
    result,output,paths=invoke(solver,tmp_path,mesh,ground,case)
    run=require_run(result,output,r'BIRTH after seepage needs an explicit initial hydraulic state for the new soil|wet BIRTH requires one fully submerged fresh soil body.*flow/Richards/structural combinations are unverified')
    solver_gates(run)
    contact_water_gates(mesh,ground,case,run)
    # At the completed contact, all owners are active and hydraulic head is
    # defined on both sides. Do not accept null/undefined or copied old values.
    state=run['frames'][-1]['state']
    pressure=np.asarray(state['pore_pressure_kpa'],dtype=float)
    assert pressure.size and np.isfinite(pressure).all()
    assert np.all(pressure>=0.)
    restart_exact(solver,tmp_path,paths,run)
