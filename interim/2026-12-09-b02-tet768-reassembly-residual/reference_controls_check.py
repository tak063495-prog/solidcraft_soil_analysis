"""Verify that diagnostics did not change the original input-only reference rows."""
import argparse
import json
from pathlib import Path
from common import read,write
from study import exact

def check(formal,controls,output):
    results=[]
    for name in ('mesh-2','mesh2-load8','mesh2-load16'):
        folder=formal/name;path=folder/('independent-all.jsonl' if (folder/'independent-all.jsonl').exists() else 'independent.jsonl')
        count=0;different=[];old_max=0.;new_max=0.
        with (controls/(name+'.jsonl')).open() as baseline,path.open() as observed:
            for a,b in zip(baseline,observed,strict=True):
                a=json.loads(a);b=json.loads(b);b.pop('linear_diagnostics',None);b.pop('source_probes',None);count+=1
                old_max=max(old_max,a['free_residual_relative']);new_max=max(new_max,b['free_residual_relative'])
                if not exact(a,b):different.append(dict(accepted_index=a['accepted_index'],fields=[k for k in a if not exact(a[k],b[k])]))
        results.append(dict(case=name,states=count,entire_original_fields_exact=not different,differing_rows=different,baseline_peak_residual=old_max,diagnostic_peak_residual=new_max,process=read(controls/(name+'.process.json'))))
    write(output,dict(comparisons=results,all_original_fields_exact=all(r['entire_original_fields_exact'] for r in results),additional_fresh_independent_histories=3,additional_accepted_states=sum(r['states'] for r in results),numerical_acceptance=False))
    assert all(r['entire_original_fields_exact'] for r in results)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('controls',type=Path);p.add_argument('output',type=Path);a=p.parse_args();check(a.formal,a.controls,a.output)
