"""Report the live-contact P/R exchange separately from the 78-path matrix."""
import argparse,json
from pathlib import Path
import numpy as np
from a02follow_binding import study,fields,metric,compare
from a02follow_run import write,sha
from a02follow_binding import HERE
from a02follow_report import load,switches,switchdiff

def report(run,main_run,main_report,out):
    out.mkdir(parents=True,exist_ok=False);paths=study.read(run/'paths.json');assert len(paths)==8 and all(p['complete'] and p['exit_code']==0 for p in paths)
    data={e:{} for e in ('native','independent')};result={}
    for p in paths:
        folder=run/p['path'];engine=p['engine'];key=f"p{p['predictor']}-r{p['repartition']}";rows=load(folder,engine);data[engine][key]=rows
        basefolder=main_run/'exchange/stress_free-cycle'/p['path']
        for name in ('mesh','ground'):assert study.read(folder/(name+'.json'))==study.read(basefolder/(name+'.json'))
        expected=study.read(basefolder/'case.json');expected['joint_system']['maximum_increment_m']=.0005
        assert expected==study.read(folder/'case.json')
        baseline=load(basefolder,engine)
        file=folder/('native/events.jsonl' if engine=='native' else 'events.jsonl');ev=[json.loads(s) for s in file.read_text().splitlines()]
        binding=study.read(folder/'binding.json');saved={Path(k).name:v for k,v in binding['source_sha256'].items()}
        for name in ('a02follow_binding.py','a02follow_run.py','a02follow_reference.py','a02follow_policies.py','a02follow_contact_exchange.py'):
            assert saved[name]==sha(HERE/name)
        if engine=='independent':force_max=study.read(folder/'independent-result.json')['maximum_force_residual']
        else:
            native_run=study.read(folder/'native/run.json');attempts=native_run['attempts']+native_run['geostatic']['attempts']
            force_max=max(a['relative_force_residual'] for a in attempts if a['accepted'])
        accepted={} 
        for r in ev:
            if (engine=='native' and r['kind']=='a02follow_iterate') or (engine=='independent' and r['kind']=='accepted'):
                accepted[r['stage'],r['progress']]=r.get('partitions',r.get('counts'))[-32:]
        assert len(accepted)==28
        counts=list(accepted.values())
        result[p['path']]=dict(work_vs_default_contact_maximum=metric([f['endpoint_work'] for _,f in rows],[f['endpoint_work'] for _,f in baseline]),
            switches_vs_default=switchdiff(switches(rows),switches(baseline)),
            maximum_reported_force_residual=force_max,
            accepted_contact_parts_minimum=min(min(x) for x in counts),accepted_contact_parts_maximum=max(max(x) for x in counts),
            accepted_contact_counts=[dict(stage=k[0],progress=k[1],counts=v) for k,v in accepted.items()],
            repartitions=sum(r['kind'] in ('repartition','a02follow_repartition') for r in ev),iterations=sum(r['kind'] in ('iterate','a02follow_iterate') for r in ev),
            predictor_special_used=sum(r.get('used',False) for r in ev if r['kind']=='predictor'))
    paired={};within={}
    for key in data['native']:
        a,b=data['native'][key],data['independent'][key]
        paired[key]=dict(fields=compare(a,b,fields.CATEGORIES),switches=switchdiff(switches(a),switches(b)))
    for e in data:
        baseline=data[e][f'p{e}-r{e}'];within[e]={}
        for k,rows in data[e].items():within[e][k]=dict(work=metric([f['endpoint_work'] for _,f in rows],[f['endpoint_work'] for _,f in baseline]),switches=switchdiff(switches(rows),switches(baseline)))
    summary=dict(passed=False,diagnostic_only=True,case='stress_free-cycle',n=4,minimum_parts=1,maximum_contact_increment_m=.0005,
        paths=paths,path_observations=result,paired=paired,within_engine=within,
        actual_contact_partition_contrast=True,original_case_geometry_loads_materials_unchanged=True,
        reused_contact_bound_is_not_a_new_physical_definition=True)
    assert all(v['accepted_contact_parts_maximum']>1 for v in result.values())
    write(out/'contact-policy-exchange.json',summary)
    final=study.read(main_report/'measurements.json')
    final['main_matrix_observer_free_compatibility_paths']=final.pop('observer_free_compatibility_paths')
    final['main_matrix_fresh_global_paths_including_compatibility']=final.pop('total_fresh_global_paths_including_compatibility')
    final['main_matrix_accepted_states']=final.pop('fresh_accepted_states')
    final['additional_actual_contact_exchange']=summary
    final['total_principal_fresh_global_paths']=len(final['paths'])+len(paths)
    final['total_observer_free_compatibility_paths']=final['main_matrix_observer_free_compatibility_paths']+sum(p['engine']=='native' for p in paths)
    final['total_fresh_global_paths_including_compatibility']=final['total_principal_fresh_global_paths']+final['total_observer_free_compatibility_paths']
    final['total_principal_accepted_states']=final['main_matrix_accepted_states']+sum(p['accepted_states'] for p in paths)
    final['total_accepted_states_including_compatibility']=final['total_principal_accepted_states']+sum(p['accepted_states'] for p in final['paths']+paths if p['engine']=='native')
    write(out/'measurements.json',final)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--main-run',type=Path,required=True);p.add_argument('--main-report',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();report(a.run,a.main_run,a.main_report,a.output)
