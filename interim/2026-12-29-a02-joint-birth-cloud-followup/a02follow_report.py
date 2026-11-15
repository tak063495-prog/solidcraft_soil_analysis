"""Fixed-order difference vectors from NEW global solves, never extrapolated."""
import argparse,json
from pathlib import Path
import numpy as np
from a02follow_binding import study,fields,metric,norm,dot,common,compare
from a02follow_run import CASES,write
QUANTITIES=('endpoint_work','plastic_dissipation','numerical_dissipation','plastic_slip','accumulated','traction','u')
def load(folder,engine):return list(fields.load_fields(folder/(engine+'-fields.jsonl')))
def values(rows,key):return np.stack([f[key] for _,f in common(rows,2)])
def switches(rows):
    result=[]
    for i,((ra,a),(rb,b)) in enumerate(zip(rows,rows[1:])):
        if not np.all(a['joint_active']) or not np.all(b['joint_active']) or not np.array_equal(a['generation'],b['generation']):continue
        for label,key in [('branch','contact_status'),('gap_open','contact_open')]:
            for gp in np.flatnonzero(a[key]!=b[key]):
                result.append(dict(kind=label,stage=rb['stage'],gp=int(gp),before=int(a[key][gp]),after=int(b[key][gp]),
                    lower=ra['progress'] if ra['stage']==rb['stage'] else 0.,upper=rb['progress'],accepted_index=i+1,generation=b['generation'].tolist()))
        for gp in np.flatnonzero(b['accumulated']>a['accumulated']):
            result.append(dict(kind='slip_activity',stage=rb['stage'],gp=int(gp),before=0,after=1,
                lower=ra['progress'] if ra['stage']==rb['stage'] else 0.,upper=rb['progress'],accepted_index=i+1,
                generation=b['generation'].tolist(),accumulated_increment_m=float(b['accumulated'][gp]-a['accumulated'][gp]),
                endpoint_status=int(b['contact_status'][gp]),exact_internal_switch_root_known=False))
    return result

def switchdiff(a,b):
    def index(rows):
        groups={}
        for r in rows:
            key=(r['kind'],r['stage'],r['gp'],r['before'],r['after']);groups.setdefault(key,[]).append(r)
        return groups
    aa,bb=index(a),index(b);matched=[];onlya=[];onlyb=[]
    for key in sorted(set(aa)|set(bb)):
        left,right=aa.get(key,[]),bb.get(key,[])
        for x,y in zip(left,right):matched.append(dict(key=list(key),left=x['upper'],right=y['upper'],difference=x['upper']-y['upper']))
        onlya.extend(left[len(right):]);onlyb.extend(right[len(left):])
    return dict(matched=len(matched),nonzero_position_differences=sum(r['difference']!=0 for r in matched),maximum_position_difference=max([abs(r['difference']) for r in matched],default=0),
        only_left=len(onlya),only_right=len(onlyb),matched_positions=matched,unmatched_left=onlya,unmatched_right=onlyb,
        convention='same stage/GP/from/to occurrence order; bracket upper accepted endpoint; no interpolation')
def report(run,out):
    out.mkdir(parents=True,exist_ok=False);paths=study.read(run/'paths.json');assert len(paths)==78 and all(p['complete'] and p['exit_code']==0 for p in paths)
    vectors={};summaries={};crosscache={};events={};exchanges={}
    for case,_,_ in CASES:
        caseout={};crosscache[case]={}
        for engine in ('native','independent'):
            grid={}
            for n in (2,4,8):
                for m in (2,4,8):
                    path=run/f'cross/n{n}-m{m}'/case/engine;grid[n,m]=load(path,engine)
                    ev=switches(grid[n,m]);events[f'cross/{case}/{engine}/n{n}-m{m}']=ev
            crosscache[case][engine]=grid
            result={}
            for q in QUANTITIES:
                y={key:values(v,q) for key,v in grid.items()};effects={}
                def effect(label,d,base):
                    vec=np.asarray(d).ravel(order='C');vectors[f'{case}/{engine}/{q}/{label}']=vec.tolist()
                    effects[label]=dict(absolute_l2=norm(vec),maximum_absolute=float(max(abs(vec),default=0)),reference_l2=norm(base),relative_l2=norm(vec)/max(norm(base),1e-12),components=len(vec))
                for n in (2,4,8):
                    for a,b in ((2,4),(4,8)):effect(f'm{a}-to-{b}_at_n{n}',y[n,b]-y[n,a],y[n,b])
                for m in (2,4,8):
                    for a,b in ((2,4),(4,8)):effect(f'n{a}-to-{b}_at_m{m}',y[b,m]-y[a,m],y[b,m])
                for na,nb in ((2,4),(4,8)):
                    for ma,mb in ((2,4),(4,8)):
                        dn=y[nb,ma]-y[na,ma];dm=y[na,mb]-y[na,ma];inter=y[nb,mb]-y[na,mb]-y[nb,ma]+y[na,ma]
                        label=f'interaction_n{na}-{nb}_m{ma}-{mb}';effect(label,inter,y[nb,mb])
                        effects[label].update(relative_to_global_effect=norm(inter)/max(norm(dn),1e-30),relative_to_local_effect=norm(inter)/max(norm(dm),1e-30),main_effect_dot=dot(dn,dm),interaction_dot_global=dot(inter,dn))
                result[q]=effects
            caseout[engine]=result
        caseout['switch_comparisons']={}
        for engine in ('native','independent'):
            grid=crosscache[case][engine];switch_rows={}
            for m in (2,4,8):
                for na,nb in ((2,4),(4,8)):
                    switch_rows[f'n{na}-to-{nb}_at_m{m}']=switchdiff(switches(grid[na,m]),switches(grid[nb,m]))
            for n in (2,4,8):
                for ma,mb in ((2,4),(4,8)):
                    switch_rows[f'm{ma}-to-{mb}_at_n{n}']=switchdiff(switches(grid[n,ma]),switches(grid[n,mb]))
            caseout['switch_comparisons'][engine]=switch_rows
        caseout['same_grid']={f'n{n}-m{m}':compare(crosscache[case]['native'][n,m],crosscache[case]['independent'][n,m],fields.CATEGORIES) for n in (2,4,8) for m in (2,4,8)}
        summaries[case]=caseout
        ex={}
        for engine in ('native','independent'):
            ex[engine]={}
            for p in ('native','independent'):
                for r in ('native','independent'):
                    key=f'p{p}-r{r}';folder=run/'exchange'/case/(engine+'-'+key);rows=load(folder,engine)
                    ex[engine][key]=rows;events[f'exchange/{case}/{engine}/{key}']=switches(rows)
        pairs={}
        for p in ('native','independent'):
            for r in ('native','independent'):
                key=f'p{p}-r{r}';left,right=ex['native'][key],ex['independent'][key]
                pairs[key]=dict(fields=compare(left,right,fields.CATEGORIES),switches=switchdiff(switches(left),switches(right)))
        effects={}
        for engine in ('native','independent'):
            baseline=ex[engine][f'p{engine}-r{engine}'];effects[engine]={}
            for key,rows in ex[engine].items():
                effects[engine][key]=dict(work=metric([f['endpoint_work'] for _,f in rows],[f['endpoint_work'] for _,f in baseline]),switches=switchdiff(switches(rows),switches(baseline)))
        baseline=metric([f['endpoint_work'] for _,f in ex['native']['pnative-rnative']],[f['endpoint_work'] for _,f in ex['independent']['pindependent-rindependent']])
        paired=[v['fields']['fields']['endpoint_work']['absolute_l2'] for v in pairs.values()]
        exchanges[case]=dict(paired=pairs,within_engine=effects,original_native_vs_independent_work=baseline,
            unexplained_paired_work_absolute_l2_min=min(paired),unexplained_paired_work_absolute_l2_max=max(paired),
            residual_is_not_a_causal_percentage=True)
    write(out/'difference-vectors.json',vectors);write(out/'branch-switches.json',events)
    # Persist actual iteration/partition observations, independently of accepted switches.
    trace={}
    for p in paths:
        folder=run/p['path'];file=folder/('native/events.jsonl' if p['engine']=='native' else 'events.jsonl')
        rows=[json.loads(s) for s in file.read_text().splitlines()]
        counts=[r.get('partitions',r.get('counts')) for r in rows if r.get('partitions',r.get('counts')) is not None]
        contact=[v[-32:] for v in counts]
        trace[p['path']]=dict(events=len(rows),repartitions=sum(r['kind'] in ('repartition','a02follow_repartition') for r in rows),
            minimum_contact_parts=min(min(x) for x in contact),maximum_contact_parts=max(max(x) for x in contact),
            iterations=sum(r['kind'] in ('iterate','a02follow_iterate') for r in rows))
    write(out/'partition-trace-summary.json',trace)
    write(out/'measurements.json',dict(passed=False,a02_complete=False,diagnostic_only=True,new_limited_acceptance=False,
        original_acceptance_all_met=False,cross=summaries,exchange=exchanges,paths=paths,
        principal_matrix_fresh_global_paths=len(paths),observer_free_compatibility_paths=sum(p['engine']=='native' for p in paths),
        total_fresh_global_paths_including_compatibility=sum(2 if p['engine']=='native' else 1 for p in paths),
        fresh_accepted_states=sum(p['accepted_states'] for p in paths),
        richardson_or_extrapolation_used=False,finest_grid_is_exact_reference=False,
        branch_switches_are_accepted_brackets_not_exact_root_locations=True,windows_ci_after_push='未確定'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();report(a.run.resolve(),a.output.resolve())
