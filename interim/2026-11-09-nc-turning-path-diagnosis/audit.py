"""Recompute diagnostics from stored stresses, loads and linear systems.

This audits native-map measurements; it does not create an independent solution.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from verify_high_order import prepare
from probe_low_preload_nc import accepted_count

def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def relative(a,b):
    a=np.asarray(a);b=np.asarray(b);d=a-b
    return dict(absolute_l2=float(np.linalg.norm(d)),relative_l2=float(np.linalg.norm(d)/max(np.linalg.norm(b),1e-30)),maximum_absolute=float(np.max(np.abs(d))))
def check_sample(sample,previous,cp,R,q,offsets,load0,holding,plane,step_size):
    states=sample['states'];u=np.array(sample['displacement_m']);u0=np.array(previous['displacement_m'])
    du=(u-u0).ravel();coords=R.T@du
    if np.max(np.abs(R@coords-du))>1e-12:raise ValueError('diagnostic constraints changed')
    if abs(plane@coords-step_size)>1e-12:raise ValueError('diagnostic displacement plane mismatch')
    force=np.zeros(len(load0));indices=[]
    for ei,eq in enumerate(q):
        ps=range(int(offsets[ei]),int(offsets[ei+1]))
        if not cp['active_elements'][ei]:
            if any(states[i]!=cp['material_states'][i] for i in ps):raise ValueError('removed material changed')
            continue
        indices.extend(ps)
        stress=np.array([states[i]['effective_stress_compression_positive_kpa'] for i in ps])
        gd=(eq['nodes'][:,None]*3+np.arange(3)).ravel()
        force[gd]-=np.einsum('qai,qa,q->i',eq['b'],stress,eq['weight'])
    target=load0-(sample['progress']-cp['stage_progress'])*holding
    r=float(np.linalg.norm(R.T@(target-force))/max(np.linalg.norm(target),1.))
    if not r<1e-8:raise ValueError('stored stress equilibrium fails original force gate')
    if abs(sample['progress']-sample['previous_progress']-sample['progress_increment'])>1e-16:
        raise ValueError('progress bookkeeping differs')
    return r,indices

def check_linear(m):
    k=np.array(m['matrix_columns']).T;b=np.array(m['rhs']);x=np.array(m['solution'])
    r=float(np.linalg.norm(k@x-b)/max(np.linalg.norm(b),1e-30))
    if not r<1e-10:raise ValueError('saved augmented linear system fails original linear gate')
    return r

def measure(evidence,output):
    output.mkdir(parents=True,exist_ok=False)
    derivative=read(evidence/'derivative-measurements.json')
    records=read(evidence/'field-records.json');seed=next(r['event'] for r in records if r['label']=='terminal-rejected:0')
    cp=seed['committed'];mesh=read(evidence/'prior/fresh-head-midpoint/input/mesh.json')
    q,offsets=prepare(mesh);nd=len(seed['reduced_residual']);R=np.zeros((len(mesh['nodes'])*3,nd))
    for g,row in enumerate(seed['dof_rows']):
        for i,c in row:R[g,i]=c
    if not np.array_equal(R.T@R,np.eye(nd)):raise ValueError('this audit requires the recorded simple DOF map')
    native=read(evidence/'prior/fresh-head-midpoint/job/recovery-checkpoint-0.json')['run']
    if native['last_committed']!=cp:raise ValueError('diagnostic seed is not the sealed saved state')
    holding=np.array(native['excavation_audits'][0]['initial_holding_load_kn']).ravel()
    load0=np.array(cp['applied_nodal_load_kn']).ravel()
    paths={};saved={}
    for name in ['path-diagnostic-coarse','path-diagnostic-fine']:
        folder=evidence/name;samples=read(folder/'samples.json');controls=read(folder/'controls.json')
        events=read(folder/'events.json');plane=np.array(controls['plane']);step=controls['step_size_m']
        if controls['force_tolerance']!=1e-8 or controls['linear_tolerance']!=1e-10:raise ValueError('diagnostic force/linear gate changed')
        if len(samples)!=64:raise ValueError('incomplete diagnostic path')
        previous=cp;residuals=[]
        for s in samples:
            r,indices=check_sample(s,previous,cp,R,q,offsets,load0,holding,plane,step)
            req=read(folder/f'call-{s["material_call"]:04d}-request.json')
            result=read(folder/f'call-{s["material_call"]:04d}-result.json')['points']
            previous_states=previous.get('states',previous.get('material_states'))
            if len(req)!=len(indices) or len(result)!=len(indices):raise ValueError('missing material call points')
            delta=np.array(s['displacement_m'])-np.array(previous['displacement_m'])
            for ip,(request,native_point) in enumerate(zip(req,result,strict=True)):
                point=indices[ip]
                if request['point']!=point or native_point['point']!=point:raise ValueError('material point ordering changed')
                if request['state']!=previous_states[point]:raise ValueError('material history origin changed')
                if native_point.get('updated_state')!=s['states'][point]:raise ValueError('material call and sample differ')
                if request['parameters']['integration']!='midpoint' or request['partitions']!=1:raise ValueError('material diagnostic controls changed')
                ei=int(np.searchsorted(offsets,point,side='right')-1);qi=point-int(offsets[ei])
                de=-q[ei]['b'][qi]@delta[q[ei]['nodes']].ravel()
                if np.max(np.abs(de-np.array(request['strain'])))>1e-13:raise ValueError('displacement and material strain differ')
            residuals.append(r);previous=s
        linear=[]
        for p in sorted(folder.glob('linear-*.json')):
            linear.append(check_linear(read(p)))
        peak=max(range(len(samples)),key=lambda i:samples[i]['progress'])
        end=samples[-1];stress=np.array([end['states'][i]['effective_stress_compression_positive_kpa'] for i in indices])
        pc=np.array([end['states'][i]['preconsolidation_pressure_kpa'] for i in indices])
        oldpc=np.array([cp['material_states'][i]['preconsolidation_pressure_kpa'] for i in indices])
        vp=np.array([sum(end['states'][i]['plastic_strain_compression_positive'][:3]) for i in indices])
        oldvp=np.array([sum(cp['material_states'][i]['plastic_strain_compression_positive'][:3]) for i in indices])
        paths[name]=dict(equilibrium_samples=len(samples),step_size_m=step,
            maximum_sampled_progress=samples[peak]['progress'],peak_sample=peak,
            end_progress=end['progress'],positive_progress_increments=sum(s['progress_increment']>0 for s in samples),
            negative_progress_increments=sum(s['progress_increment']<0 for s in samples),
            maximum_reassembled_relative_force_residual=max(residuals),
            augmented_linear_systems=len(linear),maximum_recomputed_true_linear_residual=max(linear),
            minimum_sampled_symmetric_eigenvalue_kn_m=min(s['symmetric_eigen_min_kn_m'] for s in samples),
            removed_material_frozen=True,active_material_points=len(indices),
            end_pc_decreasing_points=int(np.sum(pc<oldpc)),end_plastic_dilating_points=int(np.sum(vp<oldvp)),
            end_dry_side_points=int(np.sum(stress[:,:3].mean(1)<pc/2)),
            pc_change_range_kpa=[float((pc-oldpc).min()),float((pc-oldpc).max())],
            native_committed_increments=0,passed=False,history_converged=False)
        saved[name]=samples
    a=saved['path-diagnostic-coarse'][31];b=saved['path-diagnostic-fine'][63]
    comparison=dict(common_plane_displacement_m=.00064,progress_absolute_difference=abs(a['progress']-b['progress']),
        displacement=relative(a['displacement_m'],b['displacement_m']),fields={})
    for f in a['states'][0]:
        comparison['fields'][f]=relative([a['states'][i][f] for i in indices],[b['states'][i][f] for i in indices])
    turning=dict(coarse_sampled_peak=paths['path-diagnostic-coarse']['maximum_sampled_progress'],
        fine_sampled_peak=paths['path-diagnostic-fine']['maximum_sampled_progress'],
        peak_difference=abs(paths['path-diagnostic-coarse']['maximum_sampled_progress']-paths['path-diagnostic-fine']['maximum_sampled_progress']),
        native_next_minimum_trial=.9514002799987793,
        next_trial_above_sampled_peaks=all(.9514002799987793>p['maximum_sampled_progress'] for p in paths.values()),
        interpretation='connected native-map displacement-controlled path turns back locally; sampled peaks are not an exact critical load or proof of global nonexistence')
    contacts=[sum(p['active'] for p in r['event']['joint_report']['points']) for r in records]
    if any(contacts):raise ValueError('contact activation is outside this diagnostic scope')
    summary=dict(schema='soil-fem-nc-turning-path-diagnosis/1',passed=False,history_converged=False,
        base_commit='06bff71dcc54f7060932aff57194f04ffe0593f8',parent_commit='6429e8366ad7b02ab546cc9a9f7ff321b693962b',
        native_highest_accepted_progress=cp['stage_progress'],native_new_accepted_increments=0,
        native_accepted_increments=accepted_count(native),native_full_removal_reached=False,production_adopted=False,
        source_checkpoint_sha256=hashlib.sha256((evidence/'prior/fresh-head-midpoint/job/recovery-checkpoint-0.json').read_bytes()).hexdigest(),
        seed_zero_increment_force_checks={name:read(evidence/name/'seed-check.json') for name in paths},
        native_failure='cutback exhausted; true linear residual 1.4875587579753504e-10 (prior archived run, not rerun here)',
        native_full_solve_executed_this_task=False,native_material_replays=derivative['material_requests'],
        derivative_fields=derivative['fields'],all_derivative_fields_material_exact=True,
        contact_active_points_in_derivative_fields=contacts,paths=paths,turning=turning,
        partial_diagnostic_step_comparison=comparison,
        diagnosis='Severe MCC loading-branch tangent softening and a local fold in the connected equilibrium path. Contact switching is absent in the observed window. The final FD/implicit mismatch is small; an exact static failure load and all alternative equilibrium branches remain unproved.',
        cloud_next_steps=['bounded pseudo-arclength with native and input-only independent material evolution, to test the branch and alternatives',
            'stability and energy/dissipation audit along any alternative branch before considering a production solver change'],
        user_pc_required_for_this_small_diagnostic=False,
        specification_decision='If the intended response crosses loss of the monotone quasi-static branch, choose the physical continuation (e.g. dynamic response or a justified length/time regularization). Do not silently reverse the prescribed removal or invent parameters.',
        full_release_material_contact_history_verified=False,remaining_region_zero_release_reverified_this_task=False,
        full_release_freeze_verified=False,independent_plastic_3d_full_field_verified=False,
        full_path_step_sensitivity_verified=False,restart_to_completion_verified=False,
        windows_ci_for_new_head='unconfirmed',no_new_acceptance_thresholds=True,
        limitations=['native constitutive map, saved history and separate assembly; not an independent reference',
            'sampled local turning path is not the prescribed monotone release path',
            'no collapse/failure/stability acceptance, no global nonexistence theorem'])
    write(output/'measurements.json',summary);return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    m=measure(a.evidence.resolve(),a.output.resolve());print(json.dumps(m['turning'],indent=2))
