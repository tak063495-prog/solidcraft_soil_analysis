"""Bounded discretization/independent-history comparison for dry MCC members.

These are increment/material-integration refinements on one fixed specimen,
not a spatial mesh convergence acceptance or a certification of a CAD model.
"""
import mcc_rate_fem as rate
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np
import fem_reference as ref
from beam_reference import System as Members
from mcc_structural_reference import StructuralSystem


def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_reference_bindings(reference,source_archive=None):
    """Inputs stay exact; archived source bytes may resolve a later code edit."""
    root=Path(__file__).resolve().parents[1]
    archive=zipfile.ZipFile(source_archive) if source_archive is not None else None
    try:
        for name,digest in reference['input_and_source_sha256'].items():
            path=Path(name)
            if path.is_file() and sha(path)==digest:continue
            # An archive can resolve original repository tools only. It cannot
            # substitute a changed mesh, ground, construction or result file.
            try:relative=path.resolve().relative_to(root)
            except ValueError:relative=None
            if archive is None or relative is None or relative.parts[0]!='tools' or relative.suffix!='.py':
                raise ValueError('independent reference binding changed: '+name)
            try:data=archive.read(relative.as_posix())
            except KeyError as error:raise ValueError('original reference source is missing from archive: '+name) from error
            if hashlib.sha256(data).hexdigest()!=digest:
                raise ValueError('archived independent reference source differs: '+name)
    finally:
        if archive is not None:archive.close()
def relative(a,b,floor=1e-12):
    a,b=np.asarray(a),np.asarray(b)
    if a.shape!=b.shape:raise ValueError('independent field shape differs')
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),floor))


def physical_case(case):
    normalized=ref.normalize_stage_case(case)
    numerical={'relative_force_tolerance','relative_linear_tolerance','maximum_newton_iterations','initial_increment','minimum_increment','material_maximum_substep','tangent_perturbation','maximum_free_dofs','linear_solver','mcc_integration'}
    out={k:v for k,v in normalized.items() if k not in ('settings','verification_reject_once_stages')}
    out['physical_settings']={k:v for k,v in normalized['settings'].items() if k not in numerical}
    if out.get('beam_system') is not None:
        out['beam_system']={k:v for k,v in out['beam_system'].items() if k!='storage_budget'}
    out['stages']=[{k:v for k,v in stage.items() if k!='initial_increment'} for stage in normalized['stages']]
    return out


def material(frame,native):
    if not native:return np.asarray(frame['material_state']).reshape(-1,14)
    points=frame['state']['material_states']
    return np.array([[*p['effective_stress_compression_positive_kpa'],np.log(p['preconsolidation_pressure_kpa']),
        np.log1p(p['void_ratio']),*p['plastic_strain_compression_positive']] for p in points])


def differences(a,b,native):
    left=material(a,native);right=material(b,False)
    return dict(displacement=relative(a['state']['displacement_m'] if native else a['displacement_m'],b['displacement_m']),
        stress=relative(left[:,:6],right[:,:6]),pc=relative(np.exp(left[:,6]),np.exp(right[:,6])),
        void_ratio=relative(np.expm1(left[:,7]),np.expm1(right[:,7])),
        plastic_strain=relative(left[:,8:],right[:,8:]))


def verify(coarse_folder,fine_folder,output,source_archive=None,reference_dir=None,reference_inputs=None):
    if output.exists():raise ValueError('new comparison output required')
    reference_dir=reference_dir or coarse_folder
    paths=[coarse_folder/'ground/mesh.json',coarse_folder/'ground/ground.json',coarse_folder/'stages.json',
        coarse_folder/'stage-results.json',fine_folder/'stage-results.json',fine_folder/'ground/mesh.json',
        fine_folder/'ground/ground.json',fine_folder/'stages.json',reference_dir/'reference100.json',reference_dir/'reference200.json']
    bindings={str(p.resolve()):sha(p) for p in paths}
    if source_archive is not None:
        bindings[str(source_archive.resolve())]=sha(source_archive)
    mesh,ground,case,coarse,fine,mesh_fine,ground_fine,case_fine,r100,r200=map(read,paths)
    if mesh!=mesh_fine or ground!=ground_fine or physical_case(case)!=physical_case(case_fine):raise ValueError('refinement changed the physical problem')
    if not coarse['complete'] or not fine['complete'] or coarse['provenance']['solver_binary_sha256']!=fine['provenance']['solver_binary_sha256']:
        raise ValueError('completed runs from the same solver required')
    ref_inputs=reference_inputs or coarse_folder
    ref_paths=[ref_inputs/'ground/mesh.json',ref_inputs/'ground/ground.json',ref_inputs/'stages.json']
    for native_path,reference_path in zip(paths[:2],ref_paths[:2],strict=True):
        if sha(native_path)!=sha(reference_path):raise ValueError('mapped reference changed mesh or material/initial history')
    if physical_case(case)!=physical_case(read(ref_paths[2])):raise ValueError('mapped reference changed physical construction or kinematics')
    bindings.update({str(p.resolve()):sha(p) for p in ref_paths})
    for reference in (r100,r200):
        verify_reference_bindings(reference,source_archive)
        for name,path in zip(('mesh','ground','case'),ref_paths):
            if reference['input_and_source_sha256'][str(path.resolve())]!=sha(path):raise ValueError('reference '+name+' differs')
    if (r100['reference_steps'],r200['reference_steps'])!=(100,200):raise ValueError('100/200 reference refinement required')
    if len({len(r['frames']) for r in (coarse,fine,r100,r200)})!=1:raise ValueError('construction frame inventory differs')
    ground=ref.normalize_ground(ground);case=ref.normalize_stage_case(case)
    members=Members(mesh,case['beam_system']);projected=r200['projected_volume'];rows=[]
    for index,(a,b,c,d) in enumerate(zip(coarse['frames'],fine['frames'],r100['frames'],r200['frames'],strict=True)):
        if len({r['name'] for r in (a,b,c,d)})!=1:raise ValueError('construction stage order differs')
        for native in (a,b):
            if native['state']['active_elements']!=d['active_elements']:raise ValueError('independent active domain differs')
        row=dict(stage=a['name'],reference_refinement=differences(c,d,False),
            coarse_vs_reference=differences(a,d,True),fine_vs_reference=differences(b,d,True))
        # Fixed, explicit gates: the independent reference must refine first;
        # both Rust discretizations must remain within 1% in every field group.
        if max(row['reference_refinement'].values())>1e-4:raise ValueError('independent 100/200 path has not converged to 1e-4')
        if max(row['coarse_vs_reference'].values())>.01 or max(row['fine_vs_reference'].values())>.01:
            raise ValueError('MCC soil/member specimen differs from independent histories by more than 1%')
        kin={}
        for tag,native in (('coarse',a),('fine',b)):
            expected=d['beams'];actual=native['state']['beams']
            for key in ('active_members','active_plates'):
                if actual.get(key)!=expected.get(key):raise ValueError('independent structural activation differs')
            kin[tag]={key:relative(actual[key],expected[key],1e-9) for key in (
                'node_displacement_m','node_rotation_rad','reference_global_kinematics','plate_reference_global_kinematics')}
            if max(kin[tag].values())>.01:raise ValueError('independently propagated structural kinematics differ by more than 1%')
            response=members.solve_interface(expected,np.asarray(d['displacement_m']),native['state']['boundary_factor'])
            report=native['beam_results'];action=[];wanted=[]
            for member in report['members']:
                action.extend(np.asarray(member['local_end_actions']).reshape(-1));wanted.extend(response['members'][member['id']]['local_end_actions'])
            for plate in report.get('plates',[]):
                target=response['plates'][plate['id']]
                for sample,expected_sample in zip(plate['samples'],target['samples'],strict=True):
                    for key in ('membrane_force_kn_m','bending_moment_kn'):
                        if key not in sample:raise ValueError('unknown plate section-force format')
                        action.extend(sample[key]);wanted.extend(expected_sample[key])
            row[tag+'_member_resultant_relative_error']=relative(action,wanted,1e-6)
            if row[tag+'_member_resultant_relative_error']>.01:raise ValueError('independent member force/moment differs by more than 1%')
        row['independent_member_kinematics']=kin
        if index and case['stages'][index-1]['kind'] in ('install_beams','remove_beams','excavate'):
            stage=case['stages'][index-1];previous=r200['frames'][index-1]
            active=np.asarray(d['active_elements']);old_state=material(previous,False).reshape(len(active),-1,14)
            local=old_state[active].reshape(-1,14)
            system=StructuralSystem(mesh,ground,active,case,projected=projected)
            pressure=b['state']['surface_pressure_kpa'];factor=b['state']['boundary_factor']
            fsoil=rate.physical(mesh,ground,active,pressure,case.get('pressure_on_original_surface_only',False))[0]
            fsoil+=factor*ref.conditions_load(mesh,case['boundary'],active,skip_removed_surfaces=True)
            prepared=previous['beams']
            if stage['kind'] in ('install_beams','remove_beams'):
                prepared=members.change(prepared,np.asarray(previous['displacement_m']),stage)
            member=members.solve_interface(prepared,np.asarray(previous['displacement_m']),factor)
            if not projected and stage['kind']=='excavate':
                old_active=np.asarray(previous['active_elements']);removed=old_active & ~active
                old_pressure=fine['frames'][index-1]['state']['surface_pressure_kpa']
                old_force=rate.physical(mesh,ground,old_active,old_pressure,case.get('pressure_on_original_surface_only',False))[0]
                old_force+=factor*ref.conditions_load(mesh,case['boundary'],old_active,skip_removed_surfaces=True)
                removed_system=StructuralSystem(mesh,ground,removed,case,projected=False)
                holding=old_force.reshape(-1)-fsoil.reshape(-1)-removed_system.force(old_state[removed].reshape(-1,14))
            else:
                holding=system.force(local)-fsoil.reshape(-1)
                for node,value in member['interface_residual'].items():holding[3*node:3*node+3]+=value
            used=np.unique(ref.cells(mesh)[active]);unused=np.ones(len(mesh['nodes']),bool);unused[used]=False
            holding.reshape(-1,3)[unused]=0.
            collection=fine['excavation_audits'] if stage['kind']=='excavate' else fine['transition_audits']
            recorded=next(r for r in collection if r['stage']==stage['name'])['initial_holding_load_kn']
            row['independent_holding_relative_error']=relative(recorded,holding,1e-6)
            if row['independent_holding_relative_error']>.01:raise ValueError('independent interval holding force differs by more than 1%')
        rows.append(row)
    coarse_max={key:max(r['coarse_vs_reference'][key] for r in rows) for key in rows[0]['coarse_vs_reference']}
    fine_max={key:max(r['fine_vs_reference'][key] for r in rows) for key in coarse_max}
    for key in coarse_max:
        if coarse_max[key]>1e-8 and fine_max[key]>coarse_max[key]/2:
            raise ValueError('Rust numerical refinement did not reduce maximum '+key+' error by at least half')
    watched=paths+ref_paths+([] if source_archive is None else [source_archive])
    if bindings!={str(p.resolve()):sha(p) for p in watched}:raise ValueError('comparison inputs changed')
    report=dict(schema='soil-fem-mcc-structural-refinement/1',frames=rows,
        maximum_coarse_relative_errors=coarse_max,maximum_fine_relative_errors=fine_max,
        independent_all_history_and_member_comparison='passed',
        input_sha256=bindings,solver_sha256=coarse['provenance']['solver_binary_sha256'],
        spatial_convergence_accepted=False,arbitrary_cad_acceptance=False,goal_complete=False,p6_complete=False)
    with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in report.items() if k not in ('frames','input_sha256')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('coarse','fine','output'):p.add_argument(name,type=Path)
    p.add_argument('--source-archive',type=Path)
    p.add_argument('--references',type=Path,help='separately frozen reference directory; physical/source bindings remain mandatory')
    p.add_argument('--reference-inputs',type=Path,help='exact original independent inputs; only numerical candidate settings may differ')
    a=p.parse_args();verify(a.coarse.resolve(),a.fine.resolve(),a.output.resolve(),a.source_archive.resolve() if a.source_archive else None,
        a.references.resolve() if a.references else None,a.reference_inputs.resolve() if a.reference_inputs else None)
