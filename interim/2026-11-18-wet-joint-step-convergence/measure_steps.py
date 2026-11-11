"""All-point diagnostics on unchanged physical coordinates. No alignment/extrapolation."""
import argparse
from collections import defaultdict
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('c01_step_study',HERE/'study.py')
st=importlib.util.module_from_spec(spec);spec.loader.exec_module(st)
from study import read,write,states,STAGES,FIELDS,Water,contact
import measure as old


def entity_differences(a,b):
    ra,fa=old.load(str(a));rb,fb=old.load(str(b));ia=old.coordinates(ra);ib=old.coordinates(rb)
    keys=sorted(set(ia)&set(ib),key=lambda x:(STAGES.index(x[0]),x[1]))
    aa=[ia[k] for k in keys];bb=[ib[k] for k in keys];result={}
    for k in FIELDS:
        x=fa[k][aa];y=fb[k][bb];d=x-y
        axes=(0,)+tuple(range(2,x.ndim))
        d2=(d*d).sum(axis=axes);r2=(y*y).sum(axis=axes)
        maximum=np.max(abs(d),axis=axes)
        result[k]=[dict(index=i,difference_squared=float(v),reference_squared=float(r2[i]),maximum_absolute=float(maximum[i])) for i,v in enumerate(d2)]
    stages={}
    for s in STAGES:
        ai=[ia[k] for k in keys if k[0]==s];bi=[ib[k] for k in keys if k[0]==s]
        stages[s]={k:old.metric(fa[k][ai],fb[k][bi]) for k in FIELDS}
    return dict(common_states=len(keys),all_entities=result,stages=stages)


def attempts(run):
    return {('geostatic',a['progress_to']):a for a in run['geostatic']['attempts']} | {(a['stage'],a['progress_to']):a for a in run['attempts']}


def switches(folder,system):
    result=old.switch_data(folder,system)
    rows,_=old.load(str(folder));bykey={(r['stage'],r['progress']):i for i,r in enumerate(rows)}
    source=(read(folder.parent/'native-raw/run.json') if folder.name=='native' else None)
    am=attempts(source) if source else None
    for e in result['events']+(result['native_endpoint_status_events'] or []):
        s=e['stage'];p=e['progress_after'];r=rows[bykey[s,p]]
        i=bykey[s,p];before=rows[i-1]['factor'] if i else 0.
        e.update(boundary_factor_before=before,boundary_factor_after=r['factor'],
                 global_path_coordinate_before=STAGES.index(s)+e['progress_before'],
                 global_path_coordinate_after=STAGES.index(s)+p,
                 time_interpretation='quasistatic stage progress, not physical seconds; event bracket only')
        if am:
            a=am[s,p]
            e['attempt']=dict(accepted=a['accepted'],reason=a['reason'],newton_iterations=a['newton_iterations'],relative_force_residual=a['relative_force_residual'],maximum_true_linear_residual=max((x['true_relative_residual'] for x in a['linear_solves']),default=0.),force_limit=source['case']['settings']['relative_force_tolerance'],linear_limit=source['case']['settings']['relative_linear_tolerance'])
        else:
            e['attempt']={k:r[k] for k in ('relative_force_residual','newton_iterations','acceptance_force_tolerance','contact_partition_counts')}
    return result


def frozen_history(folder,system):
    groups=defaultdict(list)
    with (folder/'native-raw/progress.jsonl').open() as f:
        for line in f:
            v=json.loads(line)
            if v['kind']=='iterate':groups[v['stage'],v['progress']].append(v)
    assert all(len({v['committed_sha256'] for v in values})==1 for values in groups.values()),'Newton old-history mutated'
    final={key:values[-1] for key,values in groups.items()}
    maxima={'traction_kpa':0.,'plastic_slip_m':0.,'accumulated_m':0.};num=0
    previous=None;material_points=0
    from types import SimpleNamespace
    with (folder/'native-raw/accepted.jsonl').open() as f:
        for line in f:
            now=json.loads(line);key=now['stage_index'],now['stage_progress']
            assert len(now['material_states'])==324
            for m in now['material_states']:
                assert m['equivalent_plastic_strain']==m['plastic_work_kpa']==0 and m['last_return']=='elastic'
                assert not any(m['plastic_strain_compression_positive'])
            material_points+=len(now['material_states'])
            if key in final:
                assert previous is not None
                oldp=previous['joints']['points'];newp=now['joints']['points']
                counts=np.array(final[key]['partitions'][-len(oldp):])
                jump=np.array([x['jump_m'] for x in oldp]);q=np.array([x['traction_kpa'] for x in oldp]);target=np.array([x['jump_m'] for x in newp])
                qq=np.zeros_like(q);pp=np.zeros_like(q);cc=np.zeros(len(q))
                for count in np.unique(counts):
                    mask=counts==count;s=SimpleNamespace(**{k:getattr(system,k)[mask] for k in ('normal','kn','kt','mu')})
                    qi,pi,ci=contact(s,jump[mask],q[mask],target[mask],int(count))
                    qq[mask]=qi[:,0];pp[mask]=pi[:,0];cc[mask]=ci[:,0]
                maxima['traction_kpa']=max(maxima['traction_kpa'],float(np.max(abs(qq-np.array([x['traction_kpa'] for x in newp])))))
                maxima['plastic_slip_m']=max(maxima['plastic_slip_m'],float(np.max(abs(pp+np.array([x['plastic_slip_m'] for x in oldp])-np.array([x['plastic_slip_m'] for x in newp])))))
                maxima['accumulated_m']=max(maxima['accumulated_m'],float(np.max(abs(cc+np.array([x['accumulated_plastic_slip_m'] for x in oldp])-np.array([x['accumulated_plastic_slip_m'] for x in newp])))))
                num+=len(q)
            previous=now
    return dict(all_material_point_states_checked=material_points,material_equivalent_plastic_strain_and_work_zero=True,material_last_return_always_elastic=True,iterated_attempts=len(groups),all_newton_iterates_share_immutable_old_history=True,independent_replay_contact_point_increments=num,maximum_replay_differences=maxima,scope='open and recontact; final Newton per-point partition; replay uses saved inputs and is not an independent global solve')


def chord_defect(out):
    """Reassemble forces at interior points of a coarse displacement chord.

    This is a defect diagnostic, never an equilibrium solution or time-aligned field.
    The 256-step path is finite, not a continuum truth.
    """
    mesh=read(out/'mesh-h0.json');ground=read(out/'ground.json');case=read(out/'stages.json')
    op=st.prior.operators(mesh,ground,case,'cut');s=op.system
    rc,fc=old.load(str(out/'g64/fixed8'));rf,ff=old.load(str(out/'g256/fixed8'));im=old.coordinates(rf)
    zero={k:np.zeros_like(v[0]) for k,v in fc.items()}
    previous=zero;previous_fine=zero;records=[]
    for i,row in enumerate(rc):
        stage=row['stage'];p=row['progress'];start=p-1/64
        end={k:v[i] for k,v in fc.items()};fend=ff['jump'][im[stage,p]]
        old_fine_jump=previous_fine['jump'];old_fine_q=previous_fine['traction']
        curved_q=old_fine_q.copy();curved_jump=old_fine_jump.copy();curved_slip=np.zeros_like(curved_jump)
        defects=[];contact_error=[]
        for j in range(1,5):
            t=j/4;progress=start+t/64;fi=im[stage,progress]
            u=previous['u']+t*(end['u']-previous['u']);jump=np.einsum('qag,g->qa',s.c,u.ravel())
            q,flow,cum=contact(s,previous['jump'],previous['traction'],jump,2*j);q=q[:,0]
            factor=rf[fi]['factor']
            force=op.k@u.ravel()+op.seed+np.einsum('qag,qa,q->g',s.c,q,s.weight)
            load=op.base+(1-progress)*op.holding if stage=='geostatic' else op.base+factor*op.boundary
            residual=float(np.linalg.norm(op.reduce(force-load))/max(np.linalg.norm(load),1.))
            defects.append(residual)
            contact_error.append(float(np.max(abs(jump-ff['jump'][fi]))))
            qnext,pnext,_=contact(s,curved_jump,curved_q,ff['jump'][fi],8)
            curved_q=qnext[:,0];curved_slip+=pnext[:,0];curved_jump=ff['jump'][fi]
        qstraight,pstraight,_=contact(s,old_fine_jump,old_fine_q,fend,32)
        true_delta=ff['plastic_slip'][im[stage,p]]-previous_fine['plastic_slip']
        d=pstraight[:,0]-curved_slip
        records.append(dict(stage=stage,progress_before=start,progress_after=p,interior_chord_force_residuals=defects[:3],endpoint_force_residual=defects[3],maximum_chord_jump_difference_m=max(contact_error[:3]),plastic_slip_chord32_minus_curved4x8_m=d.tolist(),curved_replay_maximum_difference_from_saved_increment_m=float(np.max(abs(curved_slip-true_delta))),maximum_plastic_slip_chord_difference_m=float(np.max(abs(d)))))
        previous=end;previous_fine={k:v[im[stage,p]] for k,v in ff.items()}
    return dict(coarse_steps=64,finite_comparison_steps=256,local_partitions_chord=32,local_partitions_curved=4*8,records=records,
                stage_maxima={stage:dict(maximum_interior_force_residual=max(max(r['interior_chord_force_residuals']) for r in records if r['stage']==stage),maximum_endpoint_force_residual=max(r['endpoint_force_residual'] for r in records if r['stage']==stage),maximum_slip_chord_difference_m=max(r['maximum_plastic_slip_chord_difference_m'] for r in records if r['stage']==stage)) for stage in STAGES},
                numerical_solution=False,physical_coordinates_changed=False,continuous_reference=False)


def measure(out):
    system=Water(read(out/'mesh-h0.json'),read(out/'ground.json'),read(out/'stages.json')).contact
    pairs={};orders={};matches={};audits={};events={};entities={};frozen={}
    for kind in ('fixed8','independent','native'):
        for g in st.STEPS:
            folder=out/f'g{g}'/kind
            audits[f'{kind}/{g}']=read(folder.parent/(kind+'-audit.json'))
            events[f'{kind}/{g}']=switches(folder,system)
        for a,b in zip(st.STEPS,st.STEPS[1:]):
            key=f'{kind}/{a}__{b}';pa=out/f'g{a}'/kind;pb=out/f'g{b}'/kind
            pairs[key]=old.compare(pa,pb);entities[key]=entity_differences(pa,pb)
        for gs in zip(st.STEPS,st.STEPS[1:],st.STEPS[2:]):
            result=old.richardson([out/f'g{g}'/kind for g in gs])
            # Do not produce an extrapolated field or corrected acceptance number.
            for v in result['fields'].values():
                for k in list(v):
                    if k.startswith('richardson_correction'):del v[k]
            orders[kind+'/'+ '__'.join(map(str,gs))]=result
    for g in st.STEPS:
        matches[str(g)]=old.compare(out/f'g{g}/native',out/f'g{g}/independent')
        frozen[str(g)]=frozen_history(out/f'g{g}',system)
    interventions={}
    for name,path in [('open256-rest64',out/'g64/open256-rest64'),('open64-rest256',out/'g256/open64-rest256')]:
        for g in (64,256):interventions[f'{name}__uniform{g}']=entity_differences(path,out/f'g{g}/independent')
    result=dict(schema='soil-fem-wet-joint-global-step-diagnosis/1',starting_commit=st.BASE,passed=False,c01_decision='未達',
                global_steps=list(st.STEPS),physical_model_changed=False,tolerances_changed=False,production_solver_changed=False,
                scope='54 nodes / 12 Tet10 / 324 material GP / 32 contact GP; fixed hydrostatic CONTINUOUS / elastic MC bulk / full seven phases',
                c02_c03_c04_accepted=False,p6_complete=False,pairs=pairs,observed_orders=orders,matched_native_vs_independent=matches,path_audits=audits,
                frozen_history_replay=frozen,stage_step_interventions=interventions,
                separate_process_restart=read(out/'g256/restart-verification.json'),
                limitation='Finite 256 is not an exact reference. Nonmonotone differences, negative orders and event brackets do not create a continuous solution or a new threshold.')
    return result,events,entities,chord_defect(out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    result,events,entities,chord=measure(a.input.resolve());a.output.mkdir(parents=True,exist_ok=False)
    for name,v in [('measurements.json',result),('switch-events.json',events),('all-entity-differences.json',entities),('chord-defect.json',chord)]:write(a.output/name,v)
