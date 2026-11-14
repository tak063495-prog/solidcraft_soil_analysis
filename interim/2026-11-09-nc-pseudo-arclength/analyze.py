"""Remeasure diagnostic evidence without any material or equilibrium solve."""
import argparse,gzip,hashlib,json,re,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
from verify_high_order import prepare,physical
from low_preload_nc_reference import ClosedStickSystem
from probe_low_preload_nc import accepted_count,expected_inputs
from native_arc import read,write

def diff(a,b):
    a=np.array(a);b=np.array(b);d=a-b
    return dict(absolute_l2=float(np.linalg.norm(d)),relative_l2=float(np.linalg.norm(d)/max(np.linalg.norm(b),1e-30)),maximum_absolute=float(np.max(abs(d))))
def native_geometry(w):
    e=w/'prior';cp=read(e/'prior/fresh-head-midpoint/job/recovery-checkpoint-0.json')['run'];event=next(r['event'] for r in read(e/'field-records.json') if r['label']=='terminal-rejected:0')
    mesh=read(e/'prior/fresh-head-midpoint/input/mesh.json');q,offset=prepare(mesh);n=len(event['reduced_residual']);R=np.zeros((len(mesh['nodes'])*3,n))
    for g,row in enumerate(event['dof_rows']):
        for i,c in row:R[g,i]=c
    B=[];weight=[];pi=[]
    for ei in np.flatnonzero(event['active']):
        eq=q[ei];gd=(eq['nodes'][:,None]*3+np.arange(3)).ravel()
        for qi,v in enumerate(eq['weight']):B.append(eq['b'][qi]@R[gd]);weight.append(v);pi.append(int(offset[ei]+qi))
    holding=np.array(cp['excavation_audits'][0]['initial_holding_load_kn']);seed=cp['last_committed'];load0=np.array(seed['applied_nodal_load_kn']);lam0=seed['stage_progress']
    return SimpleNamespace(run=cp,seed=seed,R=R,B=np.array(B),weight=np.array(weight),pi=np.array(pi),holding=holding,load0=load0,lam0=lam0,mesh=mesh)
def native_force(g,states):
    stress=np.array([states[i]['effective_stress_compression_positive_kpa'] for i in g.pi]);return -np.einsum('qai,qa,q->i',g.B,stress,g.weight)
def independent_geometry(reference):
    mesh,ground,case=expected_inputs();last=read(reference/'reference.json')['last_accepted'];active=np.array(last['active_elements']);s=ClosedStickSystem(mesh,ground,active,case)
    u=read(reference/'unload.json');before=s.total_force(np.array(u['material_state'])[active].reshape(-1,14),np.array(u['displacement_m']).ravel());final=physical(mesh,ground,active,0.,True,joint_system=s.joint)[0].reshape(-1)
    return s,active,before,final-before

def compare_samples(a,b,native):
    x=a[-1];y=b[-1];out=dict(progress_absolute_difference=abs(x['progress']-y['progress']),displacement=diff(x['displacement_m'],y['displacement_m']))
    if native:
        out['material']={k:diff([p[k] for p in x['states']],[p[k] for p in y['states']]) for k in ['effective_stress_compression_positive_kpa','plastic_strain_compression_positive','preconsolidation_pressure_kpa','void_ratio']}
    else:out['material_state_14_columns']=diff(x['material_state'],y['material_state'])
    return out

def measure(w,out):
    out.mkdir(parents=True,exist_ok=False);g=native_geometry(w);run=g.run
    assert accepted_count(run)==266 and g.lam0==0.9514000415802002
    numerical=dict(native_highest_accepted_progress=g.lam0,native_accepted_increments=266,native_new_accepted_increments=0)
    origin=run['continuation']['origin'];zero_load=g.load0+g.lam0*g.holding;force=native_force(g,origin['material_states'])
    zero=dict(relative_force_residual=float(np.linalg.norm(g.R.T@zero_load-force)/max(np.linalg.norm(zero_load),1.)),
        original_reported_relative_force_residual=run['excavation_audits'][0]['zero_release_relative_force_residual'],
        active_material_points=len(g.pi),active_contact_points=0,
        removed_material_frozen_to_terminal=all(origin['material_states'][i]==g.seed['material_states'][i] for i in range(len(g.seed['material_states'])) if i not in g.pi),
        removed_contact_history_frozen_to_terminal=origin['joints']==g.seed['joints'])
    paths={};linear={};native_calls=0;native_requests=0
    for folder in sorted(w.iterdir()):
        if not folder.is_dir() or not (folder/'report.json').exists():continue
        report=read(folder/'report.json')
        if 'equilibrium_samples' in report:paths[folder.name]=report.copy()
        if (folder/'samples.json').exists():
            samples=read(folder/'samples.json');native=folder.name.startswith('native');refs=[]
            if not native:
                reference=w/('independent-input-rk4-800' if '800' in folder.name else 'independent-input-rk4' if folder.name in ('independent-arc-coarse','independent-arc-half','independent-arc-extended200') else 'independent-input-rk4-fine')
                s,active,before,delta=independent_geometry(reference)
            for row in samples:
                if native:
                    target=g.load0-(row['progress']-g.lam0)*g.holding;res=g.R.T@target-native_force(g,row['states'])
                else:
                    target=before+row['progress']*delta;res=s.reduce(target-s.total_force(np.array(row['material_state'])[active].reshape(-1,14),np.array(row['displacement_m']).ravel()))
                refs.append(float(np.linalg.norm(res)/max(np.linalg.norm(target),1.)))
            audits=[r['audit'] for r in samples]
            paths[folder.name]={**report,'maximum_reassembled_relative_force_residual':max(refs,default=None),
                'end_progress':samples[-1]['progress'] if samples else None,
                'minimum_symmetric_eigenvalue_kn_m':min((a['symmetric_eigen_min_kn_m'] for a in audits),default=None),
                'negative_directional_incremental_work_samples':sum(a['directional_incremental_work_kn_m']<0 for a in audits),
                'cumulative_internal_work_kn_m':sum(a['internal_work_kn_m'] for a in audits),
                'cumulative_external_work_kn_m':sum(a['external_work_kn_m'] for a in audits),
                'cumulative_plastic_stress_work_kn_m':sum(a['plastic_stress_work_kn_m'] for a in audits),
                'minimum_point_plastic_stress_work_kpa':min((a['minimum_point_plastic_stress_work_kpa'] for a in audits),default=None),
                'maximum_absolute_work_difference_kn_m':max((abs(a['work_difference_kn_m']) for a in audits),default=None),
                'removed_material_contact_frozen':all(a['removed_material_frozen'] and a['contact_history_frozen'] and a['active_contact_points']==0 for a in audits)}
        for name in ['linear-systems.json']:
            if not (folder/name).exists():continue
            systems=read(folder/name);residuals=[];gaps=[]
            for eq in systems:
                A=np.array(eq['matrix']);b=np.array(eq['rhs']);x=np.array(eq['solution']);r=float(np.linalg.norm(A@x-b)/max(np.linalg.norm(b),1e-30));residuals.append(r);gaps.append(abs(r-eq['true_relative_residual']))
            linear[folder.name]=dict(saved_systems=len(systems),maximum_recomputed_true_relative_residual=max(residuals,default=0.),
                maximum_reported_recomputed_difference=max(gaps,default=0.),recomputed_above_existing_1e_10=sum(r>1e-10 for r in residuals),
                note='rate-arc vectors were normalized then algebraically reconstructed before storage; this remeasurement audits the saved vectors')
        processes=folder/'processes.json'
        if folder.name.startswith('native') and processes.exists():
            for proc in read(processes):
                native_calls+=1;native_requests+=len(read(folder/Path(proc['command'][1]).name))
    def samples(n):return read(w/n/'samples.json')
    nco,nfi=samples('native-coarse'),samples('native-fine');ico,ifi=samples('independent-map-arc-coarse'),samples('independent-map-arc-fine')
    restarts={}
    for name,whole in [('native','native-coarse'),('independent','independent-arc-history400')]:
        combined=samples(name+'-restart-a')+samples(name+'-restart-b');full=samples(whole)
        restarts[name]=dict(diagnostic_only=True,whole_samples=len(full),resumed_samples=len(combined),all_saved_fields_bitwise_equal=full==combined,restart_to_full_removal_verified=False)
    branch={}
    for name in ['native-branch-search','independent-branch-search4','independent-branch-search8']:
        rows=read(w/name/'solutions.json');group=[]
        for row in rows:
            x=np.array(row['reduced_increment_m'])
            if any(np.linalg.norm(x-np.array(other['reduced_increment_m']))<1e-8 for other in group):continue
            group.append(row)
        branch[name]={**read(w/name/'report.json'),'distinct_local_root_clusters_at_1e_8_m':len(group),
            'roots':[dict(progress=r.get('target_progress',r.get('progress')),increment_norm_m=float(np.linalg.norm(r['reduced_increment_m'])),
                relative_force_residual=r['relative_force_residual'],audit=r['audit']) for r in group]}
    before_native=run['frames'];field_compare={}
    for frame in before_native:
        name=frame['name'];n=frame['state'];ind=read(w/'independent-input-rk4-800'/(name+'.json'));v=np.array(ind['material_state']).reshape(-1,14)
        field_compare[name]=dict(displacement=diff(n['displacement_m'],ind['displacement_m']),stress=diff([p['effective_stress_compression_positive_kpa'] for p in n['material_states']],v[:,:6]),
            pc=diff([p['preconsolidation_pressure_kpa'] for p in n['material_states']],np.exp(v[:,6])),void_ratio=diff([p['void_ratio'] for p in n['material_states']],np.exp(v[:,7])-1),
            plastic_strain=diff([p['plastic_strain_compression_positive'] for p in n['material_states']],v[:,8:]))
    independent_zero={}
    for name in ['independent-input-rk4','independent-input-rk4-fine','independent-input-rk4-800']:
        s,active,before,delta=independent_geometry(w/name);u=read(w/name/'unload.json');f=s.total_force(np.array(u['material_state'])[active].reshape(-1,14),np.array(u['displacement_m']).ravel())
        independent_zero[name]=dict(relative_force_residual=float(np.linalg.norm(s.reduce(f-before))/max(np.linalg.norm(before),1.)),
            note='holding independently formed from this input-only pre-cut field; lambda=0 identity, not full-release acceptance')
    original=read(w/'independent-input-rk4-800/reference.json');captured=read(w/'independent-input-rk4-800-history/reference.json')
    history=dict(read(w/'independent-input-rk4-800-history/history-binding.json'));history['terminal_and_completed_frames_exact']=original['last_accepted']==captured['last_accepted'] and all(read(w/'independent-input-rk4-800'/(f+'.json'))==read(w/'independent-input-rk4-800-history'/(f+'.json')) for f in original['completed_frames'])
    rows=0;removed_frozen=True;contact_frozen=True;old=None;first=None
    with gzip.open(w/history['history'],'rt') as stream:
        for line in stream:
            r=json.loads(line);rows+=1
            if r['stage']=='excavate_upper':
                inactive=~np.array(r['active_elements']);mat=np.array(r['material_state'])[inactive]
                if first is None:first=(mat,r['contacts'])
                removed_frozen &= np.array_equal(mat,first[0]);contact_frozen &= r['contacts']==first[1]
    history.update(recounted_rows=rows,partial_excavation_removed_material_frozen=bool(removed_frozen),partial_excavation_contact_history_frozen=contact_frozen)
    input_reports={p.name: {k:v for k,v in read(p/'reference.json').items() if k not in ['last_accepted','segments']} for p in sorted(w.glob('independent-input-*')) if p.is_dir()}
    result=dict(schema='soil-fem-nc-pseudo-arclength/1',passed=False,history_converged=False,base_commit='e21778c1e0aeeb28181322808dd0800a9747a507',parent_commit='06bff71dcc54f7060932aff57194f04ffe0593f8',
        condition=dict(nodes=54,tet10=12,material_gp=324,contact_gp=32,pc_kpa=100,preload_kpa=5,force_kn=[1,.5,0]),**numerical,
        native_full_removal_reached=False,production_adopted=False,remaining_region_lambda_zero_native=zero,remaining_region_lambda_zero_independent=independent_zero,
        paths=paths,input_only_reference_attempts=input_reports,local_root_searches=branch,linear_system_audit=linear,
        native_arc_step_comparison=compare_samples(nco,nfi,True),independent_map_arc_step_comparison=compare_samples(ico,ifi,False),
        independent_map_material_4_8_comparison=compare_samples(ico,samples('independent-map-arc-material8'),False),
        independent_local_root_4_8_comparison=diff(read(w/'independent-branch-search4/solutions.json')[-1]['material_state'],read(w/'independent-branch-search8/solutions.json')[-1]['material_state']),
        independent_vs_native_prescribed_stage_full_fields=field_compare,partial_history=history,restarts=restarts,
        native_material_processes=native_calls,native_material_point_requests=native_requests,
        peak_native_vs_independent_absolute_difference=max(r['progress'] for r in ifi)-max(r['progress'] for r in nfi),
        physical_stability='unconfirmed; loading tangent has negative directions, no dynamic or thermodynamic stability certificate',
        dissipation='plastic stress work measured; hardening free energy unspecified, total thermodynamic dissipation unconfirmed',
        alternative_branch_interpretation='two bounded local incremental roots in both maps; candidate lies on descending connected path, no disconnected/global branch theorem and no admissible monotone jump',
        full_release_material_contact_history_verified=False,full_release_freeze_verified=False,independent_plastic_3d_full_removal_verified=False,full_path_step_sensitivity_verified=False,restart_to_completion_verified=False,
        windows_ci_for_new_head='unconfirmed',no_new_physical_parameters_or_acceptance_values=True)
    write(out/'measurements.json',result);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=measure(a.evidence.resolve(),a.output.resolve());print(json.dumps(dict(passed=r['passed'],paths=len(r['paths']),peak_difference=r['peak_native_vs_independent_absolute_difference']),indent=2))
