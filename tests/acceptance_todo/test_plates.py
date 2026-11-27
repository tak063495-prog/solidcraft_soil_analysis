"""interim/2026-11-16-issue-triage/todo.ja.md B-01/B-02/B-03.

(b) = interim/2026-11-14-plate6-curved-offset/README.ja.md,
sections "オフセットの固定仕様", "曲面を受理しない理由", "検証契約と証拠".

Curved rigid-body patches are not classical shell benchmark acceptance. Valid
offset combinations get strict xfail; geometrically unspecified combinations
retain their explicit refusal assertions. No scope guard is removed.
"""
import numpy as np
import pytest
from checks import (CapabilityUnavailable, invoke, require_run, refused,
                    solver_gates, restart_exact)
from plate6_offset_study import inputs


def rigid_shell(kind, mode):
    mesh,ground,case,_=inputs(0.)
    rs=np.array([[0.,0.],[.2,0.],[0.,.2],[.1,0.],[.1,.1],[0.,.1]])
    r,s=rs.T
    points=np.c_[r,s,np.zeros(6)]
    if kind=='cylinder': points=np.c_[np.sin(r),s,np.cos(r)]
    elif kind=='sphere': points=np.c_[np.sin(r)*np.cos(s),np.sin(s),np.cos(r)*np.cos(s)]
    elif kind=='quadratic_edge': points[3,1]+=.01
    elif kind=='out_of_plane_midpoint': points[3,2]+=.01
    else: assert kind=='plane'
    cfg=case['beam_system']
    cfg['nodes']=[dict(id=str(i),position_m=p.tolist(),match_soil_node=False) for i,p in enumerate(points)]
    cfg['members']=[];cfg['nodal_loads']=[]
    cfg['translation_constraints']=[];cfg['rotation_constraints']=[]
    cfg['plate_sections'][0]['unit_weight_kn_m3']=0.
    cfg['plates']=[dict(id='rigid_shell',nodes=list(map(str,range(6))),section_id='wall_section',
                        initially_active=True,formulation='mindlin_tri6',constant_traction_kpa=[0.,0.,0.])]
    translation=np.eye(3)[mode]*.001 if mode<3 else np.zeros(3)
    rotation=np.eye(3)[mode-3]*.001 if mode>=3 else np.zeros(3)
    displacement=translation+np.cross(rotation,points)
    for i in range(6):
        for axis in range(3):
            cfg['translation_constraints'].append(dict(node_id=str(i),direction=np.eye(3)[axis].tolist(),value=float(displacement[i,axis])))
            cfg['rotation_constraints'].append(dict(node_id=str(i),direction=np.eye(3)[axis].tolist(),value=float(rotation[axis])))
    case['stages']=[];case['verification_reject_once_stages']=[]
    return mesh,ground,case,displacement,rotation


def check_rigid_shell(run,displacement,rotation):
    solver_gates(run)
    frames=run['frames'];assert frames
    for frame in frames:
        history=frame['state']['beams']
        np.testing.assert_allclose(history['node_displacement_m'],displacement,atol=1e-9,rtol=0.)
        np.testing.assert_allclose(history['node_rotation_rad'],np.tile(rotation,(6,1)),atol=1e-9,rtol=0.)
        plates=frame['beam_results']['plates'];assert len(plates)==1
        samples=plates[0]['samples'];assert samples
        for sample in samples:
            for key in ('membrane_strain','curvature_per_m','transverse_shear_strain',
                        'membrane_force_kn_m','bending_moment_kn','transverse_shear_force_kn_m'):
                np.testing.assert_allclose(sample[key],0.,atol=1e-8,rtol=0.)
        for key in ('membrane_energy_kn_m','bending_energy_kn_m','transverse_shear_energy_kn_m','drilling_penalty_energy_kn_m'):
            assert abs(plates[0][key])<=1e-8


@pytest.mark.parametrize('mode',range(6))
def test_B01_planar_rigid_patch_positive_control(solver,tmp_path,mode):
    """(b): six rigid modes; structural 1e-9 and plate raw absolute 1e-8."""
    mesh,ground,case,u,theta=rigid_shell('plane',mode)
    result,output,_=invoke(solver,tmp_path,mesh,ground,case)
    assert result.returncode==0,result.stderr
    import json
    check_rigid_shell(json.loads(output.read_text(encoding='utf-8')),u,theta)


@pytest.mark.parametrize('kind',['quadratic_edge','out_of_plane_midpoint','cylinder','sphere'])
@pytest.mark.parametrize('mode',range(6))
@pytest.mark.xfail(strict=True,raises=CapabilityUnavailable,
    reason='B-01: curved/non-straight Tri6 kernel is unimplemented; TODO B-01 / README(b), rigid patch is not shell h-convergence')
def test_B01_curved_rigid_six_modes(solver,tmp_path,kind,mode):
    """todo B-01 / (b): six rigid modes; structural 1e-9, plate raw 1e-8.

    Prescribed infinitesimal rigid motion must give zero membrane/curvature/
    shear/drilling energy using geometric normals, including curved midsides.
    """
    mesh,ground,case,u,theta=rigid_shell(kind,mode)
    result,output,_=invoke(solver,tmp_path,mesh,ground,case)
    run=require_run(result,output,r'Tri6.*(?:flat|straight|planar)|(?:flat|straight|planar).*Tri6')
    check_rigid_shell(run,u,theta)


@pytest.mark.parametrize('e',[-.05,.05])
@pytest.mark.parametrize('wet',['initial_water','staged_water'])
@pytest.mark.xfail(strict=True,raises=CapabilityUnavailable,
    reason='B-02: offset Tri6 + water is unverified; TODO B-02 / README(b), fixed plate/constraint/work gates')
def test_B02_offset_water_plate_fields_and_real_restart(solver,tmp_path,e,wet):
    """Independent neutral-plane map/recovery + restart, signed offsets.

    This criterion covers recovered plate fields and exact restart. It does
    not claim independent plastic soil/JOINT/anchor full-field acceptance.
    """
    from plate6_offset_reference import System
    from beam_reference import verify_state
    mesh,ground,case,_=inputs(e)
    if wet=='initial_water': ground['initialization']['water_level_m']=.5
    else: case['stages'].append(dict(kind='hydrostatic_water',name='wet',target_level_m=.5))
    result,output,paths=invoke(solver,tmp_path,mesh,ground,case)
    run=require_run(result,output,r'offset Tri6 with water is unverified|offset Tri6 scope excludes JOINT, anchors, water, soil BIRTH/DEATH/MATSET and unverified stage combinations')
    solver_gates(run)
    from elastic_offset import elastic_offset_gates
    elastic_offset_gates(mesh,ground,case,run)
    model=System(mesh,case['beam_system'])
    assert any(f['beam_results']['plates'] for f in run['frames'])
    for frame in run['frames']:
        state=frame['state']
        verify_state(model,state['beams'],frame['beam_results'],state['displacement_m'],state['boundary_factor'])
    restart_exact(solver,tmp_path,paths,run)


@pytest.mark.parametrize('e',[-.05,.05])
def test_B02_dry_full_field_reference_positive_control(solver,tmp_path,e):
    """todo B-02 / (b): validate the new bounded reference on accepted dry input.

    Input-only all-node u/all-GP stress relative 1e-7, structure/reference
    absolute 1e-9, plate fields/energies absolute 1e-8 and relative 1e-7.
    """
    import json
    from elastic_offset import elastic_offset_gates
    mesh,ground,case,_=inputs(e)
    result,output,_=invoke(solver,tmp_path,mesh,ground,case)
    assert result.returncode==0,result.stderr
    run=json.loads(output.read_text(encoding='utf-8'))
    solver_gates(run)
    elastic_offset_gates(mesh,ground,case,run)


@pytest.mark.parametrize('wet',['initial_water','staged_water'])
def test_B02_wet_zero_offset_legacy_limit_positive_control(solver,tmp_path,wet):
    """todo B-02 / (b): e=0 identity map validates the new wet reference.

    Solve the existing, allowed mindlin_tri6 formulation and compare the new
    reference at e=0. This does not bypass the new offset formulation's guard
    or accept nonzero offset + water. Original full-field/plate gates apply.
    """
    import copy
    import json
    from elastic_offset import elastic_offset_gates
    mesh,ground,reference_case,_=inputs(0.)
    if wet=='initial_water': ground['initialization']['water_level_m']=.5
    else: reference_case['stages'].append(dict(kind='hydrostatic_water',name='wet',target_level_m=.5))
    legacy=copy.deepcopy(reference_case)
    for plate in legacy['beam_system']['plates']:
        plate['formulation']='mindlin_tri6';plate.pop('attachment_offset_m')
    result,output,_=invoke(solver,tmp_path,mesh,ground,legacy)
    assert result.returncode==0,result.stderr
    run=json.loads(output.read_text(encoding='utf-8'))
    solver_gates(run)
    elastic_offset_gates(mesh,ground,reference_case,run)


@pytest.mark.parametrize('e',[-.0500001,.0500001])
def test_B03_outside_half_thickness_remains_refused(solver,tmp_path,e):
    """todo B-03 / (b): finite abs(e)<=t/2; do not expand the range."""
    mesh,ground,case,_=inputs(e)
    result,output,_=invoke(solver,tmp_path,mesh,ground,case)
    refused(result,output,r'offset.*(?:thickness|half)|(?:thickness|half).*offset')


@pytest.mark.parametrize('kind',['step','fold'])
def test_B03_inconsistent_shared_signed_normal_offset_remains_refused(solver,tmp_path,kind):
    """todo B-03 / (b): no averaging of inconsistent signed normal arms."""
    mesh,ground,case,_=inputs(.05)
    plate=case['beam_system']['plates'][1]
    if kind=='step': plate['attachment_offset_m']=-.05
    else: plate['nodes']=[plate['nodes'][i] for i in (0,2,1,5,4,3)]
    result,output,_=invoke(solver,tmp_path,mesh,ground,case)
    refused(result,output,r'offset.*(?:step|fold|normal)|(?:step|fold|normal).*offset')
