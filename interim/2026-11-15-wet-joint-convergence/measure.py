"""Remeasure immutable accepted fields. Diagnostic metrics never override old gates."""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import sys
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import study
from study import read, write, states, FIELDS, STAGES, Water, sha
from wet_joint_cycle_reference import compare as existing_compare


@lru_cache(maxsize=3)
def load(folder):
    rows=[]; fields={k:[] for k in FIELDS}
    for row, f in states(Path(folder)):
        rows.append(row)
        for k in fields:
            fields[k].append(f[k])
    return rows, {k:np.stack(v) for k,v in fields.items()}


def coordinates(rows):
    result={(r['stage'],r['progress']):i for i,r in enumerate(rows)}
    assert len(result)==len(rows), 'duplicate accepted coordinate'
    return result


def metric(a,b):
    assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    d=a-b; d2=float(np.sum(d*d)); r2=float(np.sum(b*b))
    return dict(difference_squared=d2,reference_squared=r2,pooled_relative_l2=float(np.sqrt(d2)/max(np.sqrt(r2),1e-12)),
                maximum_absolute=float(np.max(abs(d),initial=0)))


def compare(a,b,keys=None):
    ra,fa=load(str(a)); rb,fb=load(str(b)); ia=coordinates(ra); ib=coordinates(rb)
    selected=keys or sorted(set(ia)&set(ib),key=lambda x:(STAGES.index(x[0]),x[1]))
    assert set(selected)<=set(ia)&set(ib) and selected
    aa=np.array([ia[k] for k in selected]); bb=np.array([ib[k] for k in selected]); fields={}
    for k in FIELDS:
        x=fa[k][aa]; y=fb[k][bb]; row=metric(x,y)
        d2=((x-y).reshape(len(x),-1)**2).sum(axis=1); r2=(y.reshape(len(y),-1)**2).sum(axis=1)
        ratios=np.sqrt(d2)/np.maximum(np.sqrt(r2),1e-12); n=int(np.argmax(ratios))
        row.update(maximum_state_relative_l2=float(ratios[n]),maximum_state=dict(stage=selected[n][0],progress=selected[n][1]))
        fields[k]=row
    return dict(common_accepted_states=len(selected),finer_unmatched_states_compared=False,fields=fields,
                full_history_acceptance=False,coordinates=[dict(stage=s,progress=p) for s,p in selected])


def richardson(folders):
    data=[load(str(p)) for p in folders]; maps=[coordinates(r) for r,_ in data]
    keys=sorted(set.intersection(*(set(m) for m in maps)),key=lambda x:(STAGES.index(x[0]),x[1]))
    result={}
    for k in FIELDS:
        arrays=[f[k][[m[c] for c in keys]] for (_,f),m in zip(data,maps)]
        d0=arrays[0]-arrays[1]; d1=arrays[1]-arrays[2]; n0=float(np.linalg.norm(d0)); n1=float(np.linalg.norm(d1))
        p=float(np.log2(n0/n1)) if n0>0 and n1>0 else None
        cosine=float(np.sum(d0*d1)/(n0*n1)) if n0*n1>0 else None
        result[k]=dict(coarse_difference_l2=n0,fine_difference_l2=n1,observed_order=p,difference_direction_cosine=cosine,
                       extrapolation_is_diagnostic_only=True,accepted=False)
        if p is not None and p>0:
            correction=(arrays[2]-arrays[1])/(2**p-1)
            result[k]['richardson_correction_l2']=float(np.linalg.norm(correction))
            result[k]['richardson_correction_maximum_absolute']=float(np.max(abs(correction)))
    return dict(common_accepted_states=len(keys),fields=result,continuous_solution_certified=False,
                reason='Three levels on identical coordinates only; active-set changes and direction must be inspected. No new tolerance.')


def distribution(a,b,system):
    ra,fa=load(str(a)); rb,fb=load(str(b)); ia=coordinates(ra); ib=coordinates(rb)
    keys=sorted(set(ia)&set(ib),key=lambda x:(STAGES.index(x[0]),x[1]))
    x=fa['plastic_slip'][[ia[k] for k in keys]]; y=fb['plastic_slip'][[ib[k] for k in keys]]
    squared=(x-y)**2; total=float(squared.sum()); perpoint=squared.sum(axis=(0,2)); points=[]
    for i,v in enumerate(perpoint):
        points.append(dict(point=i,position_m=system.points[i]['position_m'],weight_m2=float(system.weight[i]),
                           difference_squared=float(v),fraction_of_difference_squared=float(v/total) if total else 0.,
                           maximum_vector_difference_m=float(np.linalg.norm(x[:,i]-y[:,i],axis=1).max()),
                           relative_l2=float(np.linalg.norm(x[:,i]-y[:,i])/max(np.linalg.norm(y[:,i]),1e-12))))
    stages={s:metric(x[[i for i,k in enumerate(keys) if k[0]==s]],y[[i for i,k in enumerate(keys) if k[0]==s]]) for s in STAGES}
    return dict(points=points,stages=stages,largest_five_points=sorted(points,key=lambda p:p['difference_squared'],reverse=True)[:5])


def switch_data(folder,system,epsilon=1e-14):
    rows,f=load(str(folder)); events=[]; previous_gap=np.zeros(len(system.points)); previous_acc=previous_gap.copy(); previous_slip=np.zeros(len(system.points),bool)
    previous_stage=None; previous_progress=0.
    for i,row in enumerate(rows):
        stage=row['stage']; progress=row['progress']
        if stage!=previous_stage:
            previous_progress=0.
        increment=f['accumulated'][i]-previous_acc
        slipping=increment>epsilon
        opened=(previous_gap<=0)&(f['gap'][i]>0)
        closed=(previous_gap>0)&(f['gap'][i]<=0)
        started=(~previous_slip)&slipping
        stopped=previous_slip&(~slipping)
        for kind,mask in (('opening',opened),('closing',closed),('slip_increment_start',started),('slip_increment_stop',stopped)):
            for point in np.flatnonzero(mask):
                e=dict(stage=stage,point=int(point),kind=kind,progress_before=previous_progress,progress_after=progress,
                       gap_before_m=float(previous_gap[point]),gap_after_m=float(f['gap'][i,point]),
                       accumulated_increment_m=float(increment[point]),factor_after=row['factor'])
                if kind in ('opening','closing'):
                    g0=previous_gap[point]; g1=f['gap'][i,point]
                    e['linear_gap_crossing_progress']=float(previous_progress-g0*(progress-previous_progress)/(g1-g0))
                events.append(e)
        previous_gap=f['gap'][i]; previous_acc=f['accumulated'][i]; previous_slip=slipping
        previous_progress=progress; previous_stage=stage
    native_events=None
    if Path(folder).name=='native':
        native_events=[]; old_status=['stick']*len(system.points); old_stage=None; old_progress=0.
        with (Path(folder).parent/'native-raw/accepted.jsonl').open(encoding='utf-8') as stream:
            for i,line in enumerate(stream):
                committed=json.loads(line);row=rows[i]
                if row['stage']!=old_stage:old_progress=0.
                now=[p['status'] for p in committed['joints']['points']]
                assert len(now)==len(old_status)
                for q,(before,after) in enumerate(zip(old_status,now)):
                    if before!=after:
                        native_events.append(dict(stage=row['stage'],point=q,from_status=before,to_status=after,
                                                  progress_before=old_progress,progress_after=row['progress'],factor_after=row['factor']))
                old_status=now;old_stage=row['stage'];old_progress=row['progress']
    normal_pressure=np.maximum(-np.einsum('tqa,qa->tq',f['traction'],system.normal),0)
    old_pressure=np.vstack((np.zeros((1,len(system.points))),normal_pressure[:-1]))
    old_gap=np.vstack((np.zeros((1,len(system.points))),f['gap'][:-1]))
    opening_mask=(old_gap<=0)&(f['gap']>0)
    drop=np.where(opening_mask,old_pressure-normal_pressure,0)
    at=np.unravel_index(int(np.argmax(drop)),drop.shape)
    coupling=dict(maximum_water_pressure_temporal_change_kpa=float(np.max(abs(f['physical_pressure']-f['physical_pressure'][0]))),
                  maximum_pore_pressure_temporal_change_kpa=float(np.max(abs(f['pore']-f['pore'][0]))),
                  maximum_effective_stress_change_between_accepted_states_kpa=float(np.max(abs(np.diff(f['stress'],axis=0)))),
                  maximum_effective_contact_pressure_drop_on_opening_interval_kpa=float(drop[at]),
                  opening_drop_location=dict(stage=rows[at[0]]['stage'],progress=rows[at[0]]['progress'],point=int(at[1])),
                  numerical_acceptance=False)
    return dict(diagnostic_slip_increment_threshold_m=epsilon,threshold_is_not_an_acceptance_tolerance=True,
                pressure_effective_stress_coupling=coupling,
                native_endpoint_status_events=native_events,
                endpoint_status_and_increment_activity_are_distinct=True,
                events=events,initial_zero_gap_is_not_an_opening_event=True,
                raw_gap_and_slip_histories_retained=True)


def aligned(a,b,ea,eb):
    """Piecewise time map anchors the first opening and slip onset per point/stage.

    No stored solution is edited; interpolation appears only in this extra metric.
    """
    ra,fa=load(str(a)); rb,fb=load(str(b)); totals=[]; refs=[]; raw=[]; eligible=0; omitted=0; shifts=[]
    for stage in STAGES:
        aa=[i for i,r in enumerate(ra) if r['stage']==stage]; bb=[i for i,r in enumerate(rb) if r['stage']==stage]
        ta=np.array([ra[i]['progress'] for i in aa]); tb=np.array([rb[i]['progress'] for i in bb])
        for q in range(fa['plastic_slip'].shape[1]):
            anchors=[]
            for kind in ('opening','slip_increment_start'):
                pair=[]
                for ev in (ea,eb):
                    rows=[x for x in ev['events'] if x['stage']==stage and x['point']==q and x['kind']==kind]
                    pair.append((rows[0]['progress_before']+rows[0]['progress_after'])/2 if rows else None)
                if all(x is not None for x in pair):
                    anchors.append(tuple(pair))
                elif pair[0]!=pair[1]:
                    omitted+=1
            anchors=sorted(set(anchors))
            xa=np.array([0.]+[x[0] for x in anchors]+[1.]); xb=np.array([0.]+[x[1] for x in anchors]+[1.])
            if np.any(np.diff(xa)<=0) or np.any(np.diff(xb)<=0):
                omitted+=1; continue
            mapped=np.interp(ta,xa,xb)
            if mapped.min()<tb.min():
                # Do not extrapolate into the missing initial state.
                keep=mapped>=tb.min()
            else:
                keep=np.ones(len(ta),bool)
            x=fa['plastic_slip'][aa,q][keep]
            y=np.column_stack([np.interp(mapped[keep],tb,fb['plastic_slip'][bb,q,c]) for c in range(3)])
            yraw=np.column_stack([np.interp(ta[keep],tb,fb['plastic_slip'][bb,q,c]) for c in range(3)])
            totals.append(x-y); refs.append(y); raw.append(x-yraw); eligible+=len(x)
            shifts.extend(float(v-u) for u,v in anchors)
    d=np.concatenate(totals); r=np.concatenate(refs); d_raw=np.concatenate(raw)
    return dict(eligible_point_states=eligible,omitted_event_or_nonmonotone_maps=omitted,
                aligned_plastic_slip=metric(d,np.zeros_like(d)) | dict(pooled_relative_l2=float(np.linalg.norm(d)/max(np.linalg.norm(r),1e-12)),reference_squared=float(np.sum(r*r))),
                unaligned_same_subset_difference_l2=float(np.linalg.norm(d_raw)),aligned_difference_l2=float(np.linalg.norm(d)),
                maximum_event_anchor_shift=float(max(map(abs,shifts),default=0.)),
                interpolation_used=True,existing_metrics_replaced=False,accepted=False)


def spatial(a,b,coarse,fine):
    ra,fa=load(str(a)); rb,fb=load(str(b)); ia=coordinates(ra); ib=coordinates(rb)
    keys=sorted(set(ia)&set(ib),key=lambda x:(STAGES.index(x[0]),x[1])); aa=[ia[k] for k in keys]; bb=[ib[k] for k in keys]
    fields={'u_common_original_nodes':metric(fa['u'][aa],fb['u'][bb,:fa['u'].shape[1]])}
    def means(f,system,children):
        parent=np.array([p['face_index']//children for p in system.points]); out=[]
        for q in range(int(parent.max())+1):
            mask=parent==q; w=system.weight[mask]; value=f[:,mask]
            out.append(np.tensordot(value,w,axes=(1,0))/w.sum())
        return np.stack(out,axis=1)
    for k in ('plastic_slip','accumulated','traction','jump','gap','physical_pressure'):
        fields[k+'_parent_face_mean']=metric(means(fa[k][aa],coarse,1),means(fb[k][bb],fine,4))
    return dict(common_accepted_states=len(keys),fields=fields,pointwise_contact_fields_compared=False,
                parent_face_mean_is_new_diagnostic_metric=True,spatial_convergence_accepted=False,
                limitation='Contact GP locations differ. Face means suppress intra-face variation; all raw GP histories remain saved. Only two h levels, no spatial order inferred.')


def measure(output):
    routes=read(output/'routes.json'); ground=read(output/'ground.json'); case=read(output/'stages.json')
    systems=[Water(read(output/f'mesh-h{h}.json'),ground,case).contact for h in (0,1)]
    failures={}; audits={}; comparisons={}; events={}; endpoint_gates={}
    for r in routes:
        route=output/r['id']
        required=['independent-result.json','independent-audit.json','independent/manifest.json']
        if r['native']:
            required+=['native-audit.json','native/manifest.json','native-raw/run.json','native-raw/verification.json','native-raw/accepted.jsonl']
        missing=[n for n in required if not (route/n).is_file()]
        if missing:
            failures[r['id']]=dict(passed=False,completed=False,missing=missing)
        if (route/'failure.json').exists():
            failures[r['id']]=read(route/'failure.json')
        for kind in ('independent','native'):
            if (route/(kind+'-audit.json')).exists():
                audits[r['id']+'/'+kind]=read(route/(kind+'-audit.json'))
                events[r['id']+'/'+kind]=switch_data(route/kind,systems[r['h']])
        if r['native'] and (route/'native-audit.json').exists():
            run=read(route/'native-raw/run.json')
            identity_path=output/f'mesh-h{r["h"]}-identity.json'
            if identity_path.exists():
                assert run['provenance']['mesh_sha256']==read(identity_path)['mesh_sha256'], 'native mesh binding differs'
            assert study.fem.normalize_ground(run['ground'])==ground
            assert study.fem.normalize_stage_case(run['case'])==study.fem.normalize_stage_case(study.configuration(case,r['g'],r['maximum']))
            verification=read(route/'native-raw/verification.json')
            assert verification['complete'] and verification['uninterrupted_run_exact_except_elapsed'] and verification['accepted_states']==7*r['g']
            comparisons[r['id']]=compare(route/'native',route/'independent')
            try:
                old_rows=existing_compare(read(route/'native-raw/run.json'),read(route/'independent-result.json'))
                endpoint_gates[r['id']]=dict(passed=True,scope='Existing displacement <0.005 and traction <0.01 at seven endpoints only',fields=old_rows)
            except AssertionError as error:
                endpoint_gates[r['id']]=dict(passed=False,reason=str(error),scope='Existing endpoint displacement/traction only')
    def p(name,kind='independent'): return output/name/kind
    global_pairs={f'{a}__{b}':compare(p(f'h0-g{a}-l8-cut'),p(f'h0-g{b}-l8-cut')) for a,b in ((16,32),(32,64),(64,128))}
    contact_pairs={f'{a}__{b}':compare(p(f'h0-g64-l{a}-cut'),p(f'h0-g64-l{b}-cut')) for a,b in ((8,16),(16,32))}
    global_fixed={f'{a}__{b}':compare(p(f'h0-g{a}-l8-cut'),p(f'h0-g{b}-l8-cut'),[(s,i/16) for s in STAGES for i in range(1,17)]) for a,b in ((16,32),(32,64),(64,128))}
    paired_global={}; paired_contact={}; spatial_rows={}; simultaneous={}; alignment={}; distributions={}
    for kind in ('independent','native'):
        for h,a,b in ((0,32,64),(0,64,128),(1,32,64)):
            x=f'h{h}-g{a}-m5e-05-cut'; y=f'h{h}-g{b}-m5e-05-cut'
            if x+'/'+kind in audits and y+'/'+kind in audits:
                paired_global[kind+f'/h{h}/{a}__{b}']=compare(p(x,kind),p(y,kind))
        for a,b in (('0.0001','5e-05'),('5e-05','2.5e-05')):
            x=f'h0-g64-m{a}-cut'; y=f'h0-g64-m{b}-cut'
            if x+'/'+kind in audits and y+'/'+kind in audits:
                paired_contact[kind+'/'+a+'__'+b]=compare(p(x,kind),p(y,kind))
        for g in (32,64):
            x=f'h0-g{g}-m5e-05-cut'; y=f'h1-g{g}-m5e-05-cut'
            if x+'/'+kind in audits and y+'/'+kind in audits:
                spatial_rows[kind+f'/g{g}']=spatial(p(x,kind),p(y,kind),*systems)
        x='h0-g32-m5e-05-cut'; y='h1-g64-m5e-05-cut'
        if x+'/'+kind in audits and y+'/'+kind in audits:
            simultaneous[kind]=spatial(p(x,kind),p(y,kind),*systems)
    for a,b in ((32,64),(64,128)):
        x=f'h0-g{a}-l8-cut'; y=f'h0-g{b}-l8-cut'; key=f'{a}__{b}'
        distributions[key]=distribution(p(x),p(y),systems[0])
        alignment[key]=aligned(p(x),p(y),events[x+'/independent'],events[y+'/independent'])
    rates=dict(global_cut=richardson([p(f'h0-g{g}-l8-cut') for g in (32,64,128)]),
               contact_cut=richardson([p(f'h0-g64-l{l}-cut') for l in (8,16,32)]))
    # All three time/contact-refined paths stay on the original mesh. This is separate from h/time refinement.
    local_simultaneous={f'{a}__{b}':compare(p(f'h0-g{a}-l{la}-cut'),p(f'h0-g{b}-l{lb}-cut')) for a,la,b,lb in ((32,4,64,8),(64,8,128,16))}
    result=dict(schema='soil-fem-wet-joint-convergence-ablation/1',starting_commit=study.BASE,passed=False,
                diagnostic_completed=not failures,failures=failures,all_history_accepted=False,p6_complete=False,
                original_16_128_parameter='global load increments per phase, not mesh subdivisions',
                mesh_levels=[dict(h=h,nodes=len(read(output/f'mesh-h{h}.json')['nodes']),tetrahedra=len(read(output/f'mesh-h{h}.json')['tetrahedra']),contact_points=len(systems[h].points)) for h in (0,1)],
                path_audits=audits,global_only_cut=global_pairs,global_fixed_16_coordinates_cut=global_fixed,
                contact_only_cut=contact_pairs,global_contact_simultaneous_cut=local_simultaneous,
                matched_native_vs_independent=comparisons,existing_endpoint_gates=endpoint_gates,matched_global_only=paired_global,
                matched_contact_cap_only=paired_contact,space_only=spatial_rows,space_time_simultaneous=simultaneous,
                volume_only=dict(independent_discrete_vs_cut=compare(p('h0-g64-m5e-05-discrete'),p('h0-g64-m5e-05-cut'))),
                richardson=rates,contact_error_distribution=distributions,event_aligned_diagnostic=alignment,
                limitations=['No all-history plastic-slip tolerance is added. Existing displacement/traction endpoint gates remain separate.',
                             'Independent cap partitions are frozen within Newton and increased after convergence if needed; native can repartition and double earlier. Same input cap is not identical partition count.',
                             'Native raw committed states are saved at every accepted checkpoint; strain, energy, reaction and pore trace in comparison are independently reassembled.',
                             'Two h levels cannot certify an asymptotic spatial sequence; 768 elements and >=256 increments are user-PC work.',
                             'Fixed hydrostatic pressure, elastic MC bulk branch and existing 7-phase full-face JOINT only.'])
    return result,events


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('evidence',type=Path); p.add_argument('output',type=Path); a=p.parse_args()
    result,events=measure(a.evidence.resolve()); a.output.mkdir(parents=True,exist_ok=False)
    write(a.output/'measurements.json',result); write(a.output/'switch-events.json',events)
