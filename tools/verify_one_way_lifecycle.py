"""Fresh-process lifecycle CLI, independent theory/FE, rollback and CAD tests."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse
import copy
from pathlib import Path
import numpy as np
from lifecycle_reference import audit,elastic
from verify_one_way_construction import read,write,invoke,restart
from verify_one_way_cases import normalized
from verify_acceptance_packet import sha
from export_one_way_construction import export


def saturated(case):
    base=case['one_way'];base['flow']['initial_head']['head_m']=1.
    for p in base['flow']['analysis']['periods']:
        p['boundaries'][0]['condition']['value']['head_m']=1.
        p['boundaries'][1]['condition']=dict(kind='outward_flux',flux_m_s=0.)
    case['transitions'][0]['action']['initial_head']['head_m']=1.


def verify(solver,output,cad=None):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False);checks={};restarts=[];elastic_checks={}
    def profile(name,mesh,ground,case,reference=False):
        folder=output/name;folder.mkdir();write(folder/'mesh.json',mesh);write(folder/'ground.json',ground);write(folder/'case.json',case)
        invoke(solver,['solve-one-way-construction',folder/'mesh.json',folder/'ground.json',folder/'case.json',folder/'result.json'],folder/'solve.log')
        run=read(folder/'result.json');checks[name]=audit(mesh,run);write(folder/'verification.json',checks[name]);export(folder/'mesh.json',folder/'result.json',solver,folder/'results')
        if reference:elastic_checks[name]=elastic(mesh,run);write(folder/'elastic-comparison.json',elastic_checks[name])
        return run,folder
    ground=read(root/'examples/ground-one-way-lifecycle.json');case=read(root/'examples/one-way-lifecycle.json');mesh=read(root/'examples/seepage-column-tet4-mesh.json')
    for order in ('tet4','tet10'):
        m=read(root/f'examples/seepage-column-{order}-mesh.json')
        for policy in ('preserve_head','preserve_water'):
            c=copy.deepcopy(case);c['transitions'][1]['action']['water']=dict(kind=policy)
            if policy=='preserve_water':c['transitions'][1]['action']['water']['reference_head']=dict(head_m=0.)
            run,folder=profile(f'{order}-{policy}',m,ground,c,True);restarts+=restart(solver,folder,run,[1,2,3,4,5,6])
    # Closed-form saturated, linear elastic one-dimensional birth on both orders.
    theory=[]
    for order in ('tet4','tet10'):
        m=read(root/f'examples/seepage-column-{order}-mesh.json');c=copy.deepcopy(case);saturated(c);c['transitions']=c['transitions'][:1];c['one_way']['flow']['analysis']['periods']=c['one_way']['flow']['analysis']['periods'][:2];c['one_way']['exterior_periods']=c['one_way']['exterior_periods'][:2]
        run,folder=profile(f'{order}-saturated',m,ground,c,True);expected=-(48.75+.375*(16.+10.*.7/1.7))/12000.;actual=np.mean(np.asarray(run['last_committed']['mechanical']['displacement_m'])[m['node_sets']['top'],2]);assert abs(actual-expected)<2e-10
        assert abs(run['transitions'][0]['water']['water_added_m3']-.5*.7/1.7)<1e-13
        theory.append(dict(element=order,top_displacement_m=float(actual),closed_form_m=expected,absolute_error_m=abs(float(actual)-expected)))
    c=copy.deepcopy(case);saturated(c);c['transitions']=c['transitions'][:1];c['one_way']['flow']['analysis']['periods']=c['one_way']['flow']['analysis']['periods'][:2];c['one_way']['exterior_periods']=c['one_way']['exterior_periods'][:2]
    for m in c['one_way']['flow']['materials']:m['retention']['saturated_specific_storage_per_m']=.1
    run,_=profile('positive-storage',mesh,ground,c,True);assert abs(run['transitions'][0]['water']['water_added_m3']-(.5*.7/1.7+.1*.125))<1e-13
    # A genuinely plastic candidate is rejected after material history changes.
    plastic=copy.deepcopy(ground)
    for m in plastic['materials'][1:]:m['mohr_coulomb']['cohesion_kpa']=2.;m['mohr_coulomb']['friction_deg']=10.
    c=copy.deepcopy(case);c['verification_reject_once_transitions']=['place','improve'];run,folder=profile('plastic-cutback',mesh,plastic,c);restarts+=restart(solver,folder,run,[2,3,5])
    assert any(not a['accepted'] and a['plastic_candidate_changed'] for a in run['attempts'])
    assert all(a['committed_pair_sha256_before']==a['committed_pair_sha256_after'] for a in run['attempts'] if not a['accepted'])
    clean=copy.deepcopy(c);clean['verification_reject_once_transitions']=[]
    for t in clean['transitions']:t['initial_increment']=.25
    direct,_=profile('plastic-clean',mesh,plastic,clean);assert [f['state'] for f in run['frames']]==[f['state'] for f in direct['frames']]
    # Time-zero MCC cut -> inactive MATSET -> fresh MC BIRTH retains MCC archive.
    mcc=read(root/'examples/ground-one-way-mcc.json')
    for layer in mcc['layers']:layer['consolidation']=dict(definition='mcc_pc_over_mean_effective',ratio=2.)
    mcc['materials'].append(copy.deepcopy(ground['materials'][2]));mcc['construction_material_ids']=[3]
    c=copy.deepcopy(case);c['initially_inactive_bodies']=[]
    for p in c['one_way']['flow']['materials'][:2]:p['retention']['saturated_water_content']=.6
    c['one_way']['flow']['materials'][2]['retention']['saturated_water_content']=.7/1.7
    c['excavations']=[dict(name='cut',before_period='hold',bodies=['UpperClay'],water_boundary_policy='explicit_period_boundaries',cut_traction_policy='explicit_period_reservoirs')]
    born=copy.deepcopy(c['transitions'][0]);born['before_period']='changed';t=copy.deepcopy(c['transitions'][1]);t['name']='replace_inactive';t['before_period']='born';t['action']['history']=dict(mode='replace',states=born['action']['states']);t['action']['water']=dict(kind='inactive');c['transitions']=[t,born];c['one_way']['flow']['analysis']['periods'][1]['boundaries'][1]['face_set']='floor'
    run,folder=profile('mcc-cut-replace-birth',mesh,mcc,c);restarts+=restart(solver,folder,run,[1,2,3,5,7,8])
    assert run['excavations'][0]['time_s']==0. and run['excavations'][0]['removed_material_points']
    assert all('preconsolidation_pressure_kpa' in p['state'] for p in run['excavations'][0]['removed_material_points'])
    assert all('preconsolidation_pressure_kpa' not in run['last_committed']['mechanical']['material_states'][p['integration_point_index']] for p in run['excavations'][0]['removed_material_points'])
    # Omitted derived points reconstruct using changed material water laws.
    c=copy.deepcopy(case);c['one_way']['store_transfer_integration_points']=False;run,_=profile('omitted',mesh,ground,c)
    ordinary=read(output/'tet4-preserve_water/result.json')
    def primary(f):
        out=copy.deepcopy(f['state']);out['fluid']['hydraulic_points']=[];out['fluid']['mechanical_points']=[];return out
    assert [primary(f) for f in ordinary['frames']]==[primary(f) for f in run['frames']]
    # Explicit future coverage is checked before creating any result.
    c=copy.deepcopy(case);c['one_way']['flow']['analysis']['periods'][1]['boundaries'].pop();write(output/'invalid-case.json',c)
    invoke(solver,['solve-one-way-construction',output/'omitted/mesh.json',output/'omitted/ground.json',output/'invalid-case.json',output/'invalid-result.json'],output/'invalid.log',False);assert not (output/'invalid-result.json').exists()
    # A genuine incompatible water transfer fails at its actual committed state.
    c=copy.deepcopy(case);c['one_way']['flow']['materials'][2]['retention']['residual_water_content']=.45;write(output/'failure-case.json',c)
    invoke(solver,['solve-one-way-construction',output/'omitted/mesh.json',output/'omitted/ground.json',output/'failure-case.json',output/'failure-result.json'],output/'failure.log',False)
    failed=read(output/'failure-result.json');assert not failed['complete'] and failed['last_committed']==failed['frames'][-1]['state'] and failed['last_committed']['hydraulic']['time_s']==200. and len(failed['transitions'])==1
    assert failed['attempts'][-1]['committed_pair_sha256_before']==failed['attempts'][-1]['committed_pair_sha256_after']
    # The actual SolidCraft-generated mesh has a named conformal layer interface.
    if cad:
        m=read(cad/'mesh.json');g=read(root/'examples/ground-local.json');c=copy.deepcopy(case)
        for mat in g['materials']:
            mat.pop('mcc',None);mat['unit_weight_kn_m3']=16.;mat['mohr_coulomb']=dict(young_kpa=10000.,poisson=.25,cohesion_kpa=500.,friction_deg=30.,dilation_deg=0.,tension_cutoff_kpa=100.,flow='non_associated')
        for layer in g['layers']:layer['consolidation']=dict(definition='not_applicable');layer['void_ratio']=dict(mode='explicit',value=.7)
        g['materials'].append(copy.deepcopy(ground['materials'][2]));g['construction_material_ids']=[3]
        c['one_way']['mechanics']['settings']['kinematics']='side_rollers';c['one_way']['mechanics']['pressure_on_original_surface_only']=True;c['one_way']['mechanics']['monitors']={};c['one_way']['mechanics']['monitor_sets']={}
        for i,p in enumerate(c['one_way']['flow']['analysis']['periods']):p['boundaries'][1]['face_set']='interface:LowerClay:UpperClay' if i==0 else 'free_surface'
        folder=output/'cad-wrapper';folder.mkdir();write(folder/'ground.json',g);write(folder/'case.json',c)
        import subprocess
        args=['pwsh','-NoProfile','-File',str(root/'scripts/Invoke-OneWayConstruction.ps1'),'-Design',str(root/'examples/local-excavation-cylinder.json'),'-ElementOrder','2','-MeshSizeM','4','-GroundFile',str(folder/'ground.json'),'-CaseFile',str(folder/'case.json'),'-SolverPath',str(solver),'-OutputDirectory',str(folder/'run')]
        from verify_checkpoint import invoke as execute
        execute(args,folder/'wrapper.log',0);native=read(folder/'run/result.json');captured=read(folder/'run/mesh.json')
        def geometry_only(value):
            out=copy.deepcopy(value)
            for r in out['regions']:r.pop('source_sha256')
            out.pop('nodes')
            return out
        assert geometry_only(captured)==geometry_only(m)
        coordinate_difference=float(np.max(np.abs(np.asarray(captured['nodes'])-m['nodes'])))
        assert coordinate_difference<1e-11
        checks['cad-wrapper']=audit(captured,native);elastic_checks['cad-wrapper']=elastic(captured,native);write(folder/'verification.json',checks['cad-wrapper']);write(folder/'elastic-comparison.json',elastic_checks['cad-wrapper'])
        write(folder/'cad-input-link.json',dict(fresh_cli_capture=True,topology_sets_and_materials_equal_to_parent_capture=True,maximum_coordinate_difference_m=coordinate_difference,coordinate_comparison_tolerance_m=1e-11,source_mesh_sha256=sha(cad/'mesh.json'),analysis_mesh_sha256=sha(folder/'run/mesh.json')))
    result=dict(schema='soil-fem-one-way-lifecycle-cli-verification/1',passed=True,solver_sha256=sha(solver),cad_exercised=bool(cad),restarts=restarts,profiles=checks,elastic_comparisons=elastic_checks,saturated_birth_theory=theory,plastic_rollback_clean_path_exact=True,mcc_archive_survives_inactive_replacement_and_birth=True,invalid_future_boundary_rejected=True,genuine_infeasible_water_transfer_preserves_last_committed=True,omitted_primary_exact=True)
    write(output/'verification.json',result);print('Lifecycle CLI verified:',len(checks),'profiles,',len(restarts),'separate-process restarts')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',required=True,type=Path);p.add_argument('--output-dir',required=True,type=Path);p.add_argument('--cad',type=Path);a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve(),a.cad.resolve() if a.cad else None)
