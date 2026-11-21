"""Fixed-order fresh matrix comparison and exhaustive category boundary records."""
import argparse
from collections import Counter
import itertools
from c04follow_common import *

FIELDS=('plastic_dissipation','traction','plastic_slip','u')

def metric(a,b):
    result=compare(a,b);result['denominator']=max(result['reference_l2'],1e-12)
    result['denominator_floor']=1e-12
    return result

def dot(a,b):
    return math.fsum(float(x)*float(y) for x,y in zip(a,b,strict=True))

def point_observation(fields,previous,gp,normal,decision):
    q=fields['traction'][gp];qn=dot(q,normal)
    tau=[float(x)-qn*float(n) for x,n in zip(q,normal,strict=True)]
    radius=float(fields['friction_limit'][gp]);shear=norm(tau)
    increment=[float(x)-float(y) for x,y in zip(fields['plastic_slip'][gp],previous[gp],strict=True)]
    margin=radius-shear;scaled=margin/(1.+radius+shear)
    return dict(status=int(fields['contact_status'][gp]),returned_radius_kpa=radius,returned_shear_norm_kpa=shear,
        returned_margin_kpa=margin,returned_margin_scaled=scaled,returned_in_observation_band=abs(scaled)<=1e-14,
        saved_plastic_increment_m=increment,saved_plastic_increment_exact_zero=all(x==0. for x in increment),
        actual_classification_operands=decision)

def classify(a,b):
    both_band=a['returned_in_observation_band'] and b['returned_in_observation_band']
    both_zero=a['saved_plastic_increment_exact_zero'] and b['saved_plastic_increment_exact_zero']
    if both_band and both_zero:return 'returned_boundary_and_both_saved_increments_exact_zero'
    if both_band:return 'returned_boundary_with_nonzero_saved_increment; inspect_raw_magnitude'
    return 'outside_returned_roundoff_band; substantive_state_difference_not_ruled_out'

def load_path(folder,m,engine):
    data=list(rows(folder/'fields.jsonl'));n=7*m+192
    if len(data)!=n:raise ValueError(('accepted count',folder,len(data),n))
    case=read(folder/'case.json')
    phases=[dict(name='geostatic',n=m)]+[dict(name=s['name'],n=round(1/s['initial_increment'])) for s in case['stages']]
    coords=[(s['name'],i/s['n']) for s in phases for i in range(1,s['n']+1)]
    if [(v['row']['stage'],v['row']['progress']) for v in data]!=coords:raise ValueError('accepted schedule')
    if engine=='native':
        raw=list(rows(folder/'raw/accepted.jsonl'));diags=list(rows(folder/'raw/constitutive.jsonl'))
        parts=[v['partitions'][-6:] for v in diags]
        decisions=[];last_decisions=[None]*6;old_jumps=[[0.,0.,0.] for _ in range(6)]
        for index,state in enumerate(raw):
            current=[]
            for gp,point in enumerate(state['joints']['points']):
                d=dict(point['c04follow_decision'])
                if d['kind']=='unchanged':
                    d['shortcut_condition']=dict(lhs=point['jump_m'],rhs=old_jumps[gp],operator='exact array ==',result=point['jump_m']==old_jumps[gp])
                    if not d['shortcut_condition']['result']:raise ValueError('native unchanged shortcut operands differ')
                    d['inherited_decision']=last_decisions[gp]
                else:
                    last_decisions[gp]=dict(index=index,stage=data[index]['row']['stage'],progress=data[index]['row']['progress'],decision=d.copy())
                current.append(d);old_jumps[gp]=point['jump_m']
            decisions.append(current)
        result=read(folder/'raw/run.json');assert result['complete']
        from joint_reference import System
        system=System(read(folder.parent/'mesh.json'),case['joint_system'],result['joint_integration_points'])
        layout=result['joint_integration_points']
        maximum=max(v['relative_force_residual'] for v in diags)
    else:
        result=read(folder/'result.json');assert result['complete'] and result['reference_uses_native_results'] is False
        parts=[v['contact_partitions'] for v in result['rows']]
        decisions=[v['c04follow_decisions'] for v in result['rows']]
        from joint_reference import System
        system=System(read(folder.parent/'mesh.json'),case['joint_system'])
        layout=system.points
        maximum=max(v['row']['relative_force_residual'] for v in data)
    if len(parts)!=n or any(len(v)!=6 for v in parts):raise ValueError('partition count')
    return dict(data=data,vectors={k:np.asarray([v['fields'][k] for v in data]) for k in FIELDS},
        parts=parts,decisions=decisions,normal=system.normal,layout=layout,maximum_force_residual=maximum,
        original_force_tolerance=case['settings']['relative_force_tolerance'] if engine=='native' else 1e-10)

def pair(a,b,output,name):
    if [(v['row']['stage'],v['row']['progress']) for v in a['data']]!=[(v['row']['stage'],v['row']['progress']) for v in b['data']]:raise ValueError('pair schedule')
    if not np.array_equal(a['normal'],b['normal']):raise ValueError('normal order')
    previous_a=np.zeros((6,3));previous_b=previous_a.copy();mismatches=[]
    for index,(aa,bb) in enumerate(zip(a['data'],b['data'],strict=True)):
        fa,fb=aa['fields'],bb['fields']
        for gp in range(6):
            if fa['contact_status'][gp]==fb['contact_status'][gp]:continue
            left=point_observation(fa,previous_a,gp,a['normal'][gp],a['decisions'][index][gp])
            right=point_observation(fb,previous_b,gp,b['normal'][gp],b['decisions'][index][gp])
            mismatches.append(dict(index=index,stage=aa['row']['stage'],progress=aa['row']['progress'],gp=gp,
                native_parts=a['parts'][index][gp],independent_parts=b['parts'][index][gp],
                native=left,independent=right,classification=classify(left,right),
                category_mechanism=('native exact-zero-jump shortcut retains preceding status; independent reclassifies current last_flow'
                    if left['actual_classification_operands'].get('kind')=='unchanged' else 'inspect actual operands')) )
        previous_a=np.asarray(fa['plastic_slip']);previous_b=np.asarray(fb['plastic_slip'])
    write(output/(name+'.json'),mismatches)
    return dict(fields={k:metric(a['vectors'][k],b['vectors'][k]) for k in FIELDS},
        category_mismatches=len(mismatches),classifications=dict(Counter(v['classification'] for v in mismatches)),
        category_mechanisms=dict(Counter(v['category_mechanism'] for v in mismatches)),
        both_returned_in_band=sum(v['native']['returned_in_observation_band'] and v['independent']['returned_in_observation_band'] for v in mismatches),
        both_saved_increments_exact_zero=sum(v['native']['saved_plastic_increment_exact_zero'] and v['independent']['saved_plastic_increment_exact_zero'] for v in mismatches),
        maximum_scaled_return_margin=max((abs(v[e]['returned_margin_scaled']) for v in mismatches for e in ('native','independent')),default=0.),
        maximum_saved_increment_m=max((abs(x) for v in mismatches for e in ('native','independent') for x in v[e]['saved_plastic_increment_m']),default=0.),
        exhaustive_mismatch_file=name+'.json')

def factorial(paths,output,name):
    # Differences are I minus N. Means use fsum in a fixed two-term order.
    result={};signed={}
    for field in FIELDS:
        nn,ni,inn,ii=[paths[p+r]['vectors'][field] for p,r in (('N','N'),('N','I'),('I','N'),('I','I'))]
        pn=inn-nn;pi=ii-ni;rn=ni-nn;ri=ii-inn
        pmean=np.asarray([math.fsum((float(x),float(y)))*.5 for x,y in zip(pn.ravel(),pi.ravel(),strict=True)]).reshape(nn.shape)
        rmean=np.asarray([math.fsum((float(x),float(y)))*.5 for x,y in zip(rn.ravel(),ri.ravel(),strict=True)]).reshape(nn.shape)
        interaction=pi-pn
        values=dict(P_at_RN=pn,P_at_RI=pi,R_at_PN=rn,R_at_PI=ri,P_mean=pmean,R_mean=rmean,PR_mixed_difference=interaction)
        signed[field]={k:v for k,v in values.items()}
        result[field]={k:metric(v,np.zeros_like(v))|dict(reference_for_effect_l2=norm(nn),relative_to_NN=norm(v)/max(norm(nn),1e-12)) for k,v in values.items()}
    write(output/(name+'.json'),signed)
    return dict(fields=result,signed_file=name+'.json',difference_sign='I minus N; PR = II-IN-NI+NN, implemented (II-NI)-(IN-NN)',norms_are_not_additive_causal_percentages=True)

def report(run,output):
    sys.path.insert(0,str(ROOT/'tools'))
    output.mkdir(parents=True,exist_ok=False)
    summary=dict(passed=False,c04_complete=False,diagnostic_only=True,reduction_order='accepted coordinate, fixed GP/node, component C order; math.fsum',
        observation_band=1e-14,observation_band_is_acceptance=False,levels={},fresh_equilibrium_paths=16,fresh_accepted_states=7552,
        mechanical32='A-04 saved values quoted separately only; not loaded or reaggregated here',rates_calculated=False,extrapolation=False)
    for m in (16,64):
        paths={e:{p+r:load_path(run/f'm{m}-{e}-p{p}r{r}',m,e) for p,r in itertools.product('NI',repeat=2)} for e in ('native','independent')}
        pairs={}
        for lp,lr,rp,rr in itertools.product('NI',repeat=4):
            left=lp+lr;right=rp+rr;name=f'm{m}-native{left}-independent{right}'
            pairs[left+'-'+right]=pair(paths['native'][left],paths['independent'][right],output,name)
        effects={e:factorial(paths[e],output,f'm{m}-{e}-signed-effects') for e in paths}
        counts={e:{k:dict(sorted(Counter(str(x) for row in v['parts'] for x in row).items())) for k,v in paths[e].items()} for e in paths}
        audit={e:{k:dict(accepted_states=len(v['data']),maximum_force_residual=v['maximum_force_residual'],original_force_tolerance=v['original_force_tolerance'],all_contact_updates=len(v['parts'])*6) for k,v in paths[e].items()} for e in paths}
        write(output/f'm{m}-gp-layout.json',{e:paths[e]['NN']['layout'] for e in paths})
        summary['levels'][str(m)]=dict(original=pairs['NN-II'],fully_exchanged=pairs['II-NN'],same_policy={p+r:pairs[p+r+'-'+p+r] for p,r in itertools.product('NI',repeat=2)},all_cross_policy_pairs=pairs,effects=effects,contact_partition_counts=counts,solver_audits=audit)
    write(output/'measurements.json',summary)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();report(a.run.resolve(),a.output.resolve())
