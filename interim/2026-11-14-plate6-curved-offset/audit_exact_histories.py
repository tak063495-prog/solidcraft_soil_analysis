"""Sign-of-zero-aware exact audit of the recorder's complete JSON histories."""
import argparse
import copy
import json
from pathlib import Path


def without_clock(v):
    if isinstance(v,dict):return {k:without_clock(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [without_clock(x) for x in v]
    return v


def canonical(v):
    return json.dumps(without_clock(v),sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False)


def accepted(v):
    v=copy.deepcopy(v)
    v['attempts']=[a for a in v['attempts'] if a['accepted']]
    v['geostatic']['attempts']=[a for a in v['geostatic']['attempts'] if a['accepted']]
    v['geostatic'].pop('verification_rejection_done',None)
    v['case'].pop('verification_reject_once_stages',None)
    v['case']['settings'].pop('initial_increment')
    for stage in v['case']['stages']:stage.pop('initial_increment',None)
    v['ground']['initialization']['geostatic'].pop('verification_reject_once',None)
    for key in ('ground_sha256','case_sha256'):v['provenance'].pop(key)
    return v


def audit(folder):
    read=lambda p:json.loads(p.read_text(encoding='utf-8'))
    rows=[]
    for name in ('positive','negative'):
        p=folder/name;run=read(p/'stage-results.json');assert run['complete']
        clean=read(p/'no-rejection-run.json');assert clean['complete']
        assert canonical(accepted(run))==canonical(accepted(clean)),(name,'accepted Run')
        resumes=sorted(p.glob('resumed-*.json'));assert len(resumes)==5
        for r in resumes:assert canonical(run)==canonical(read(r)),str(r)
        for cache in (0,64):
            bounded=read(p/f'bounded-{cache}-run.json');assert bounded['complete']
            for key in ('frames','last_committed','attempts','geostatic','transition_audits'):
                assert canonical(run[key])==canonical(bounded[key]),(name,cache,key)
        rows.append(dict(case=name,accepted_run_exact=True,restarts_exact=len(resumes),cache_variants_exact=2))
    a,b=[read(folder/'zero'/name) for name in ('offset.json','planar.json')]
    assert a['complete'] and b['complete']
    for key in ('frames','last_committed','attempts','geostatic','transition_audits'):
        assert canonical(a[key])==canonical(b[key]),('zero',key)
    return dict(passed=True,signed_zero_preserved=True,restarts_exclude_only=['elapsed_seconds'],rows=rows,
                zero_offset_physical_history_exact=True,rejection_comparison='All accepted Run data, separating only rejected attempts, rejection-control flags/initial-step inputs and their input hashes; elapsed_seconds excluded.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args();result=audit(args.folder)
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,ensure_ascii=False);f.write('\n')
    print(json.dumps(result,indent=2))
