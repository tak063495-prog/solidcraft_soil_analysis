"""All original residuals plus explicit high-precision observation selections."""
import argparse
import json
import numpy as np
from pathlib import Path
from common import read, write
from precision import Probe
from field_probe import fields

def analyze(formal,out):
    out.mkdir(parents=True,exist_ok=False)
    plan=read(formal/'fresh-solves.json')['cases']
    results=[];field_results={}
    coarse=(read(formal/'beam-positive/mesh.json'),read(formal/'beam-positive/stages.json'))
    for case in plan:
        folder=formal/case['name'];comparison=read(folder/'comparison.json')['all_point_component_comparisons']
        details={}
        for kind in ('native','independent'):
            metric=kind+'_reassembled_relative_residual'
            maximum=max(comparison,key=lambda r:r[metric])['accepted_index']
            selected={maximum}
            if case['level']==2 and case['steps']==4:selected.update(r['accepted_index'] for r in comparison)
            # Always inspect the loaded endpoint as well as each original maximum.
            probe=Probe(folder,kind=='native');rows=[];linear=[]
            source=folder/('native/accepted.jsonl' if kind=='native' else 'independent.jsonl')
            with source.open() as f:
                for line in f:
                    row=json.loads(line)
                    if 'linear_diagnostics' in row:linear.append(dict(accepted_index=row['accepted_index'],stage_index=row['stage_index'],progress=row['progress'],**row['linear_diagnostics']))
                    if row['accepted_index'] in selected or (row['stage_index']==2 and row['progress']==1.):rows.append(probe.observe(row))
            details[kind]=dict(selected_accepted_indices=[r['accepted_index'] for r in rows],original_maximum_accepted_index=maximum,observations=rows,all_direct_linear_observations=linear)
            print(case['name'],kind,'precision',len(rows),flush=True)
        run=read(folder/'native/run.json')
        details['native_solver_attempts']=dict(settings=read(folder/'stages.json')['settings'],geostatic=run.get('geostatic'),attempts=run['attempts'])
        # Keep full Run in raw evidence; only actual attempts/settings in this summary.
        if details['native_solver_attempts']['geostatic'] is not None:
            geo=details['native_solver_attempts']['geostatic'];details['native_solver_attempts']['geostatic']={k:v for k,v in geo.items() if k in ('attempts','progress','settings')}
        write(out/(case['name']+'-residual.json'),dict(passed=False,case=case['name'],details=details))
        result=fields(folder,coarse);write(out/(case['name']+'-fields.json'),result);field_results[(case['level'],case['steps'])]=result
        results.append(dict(name=case['name'],level=case['level'],steps=case['steps'],residual_file=case['name']+'-residual.json',fields_file=case['name']+'-fields.json'))
    pairs=[]
    def delta(a,b):
        assert [(r['stage_index'],r['progress']) for r in a]==[(r['stage_index'],r['progress']) for r in b]
        sums={};points=[]
        for aa,bb in zip(a,b,strict=True):
            row={}
            for key in ('shear_recovered_global_kn_m','shear_direct_global_kn_m','shear_unprojected_global_kn_m','work_sum_kn_m','nodal_source_shear_global_kn_m','other_source_shear_global_kn_m'):
                if key not in aa:continue
                d=np.asarray(aa[key])-bb[key];row[key]=float(np.linalg.norm(d));sums[key]=sums.get(key,0.)+float(np.sum(d*d))
            work={key:aa['work_components_kn_m'][key]-bb['work_components_kn_m'][key] for key in aa['work_components_kn_m']}
            points.append(dict(stage_index=aa['stage_index'],progress=aa['progress'],differences=row,work_component_signed_differences_kn_m=work))
        return dict(absolute_l2_all_points={k:float(np.sqrt(v)) for k,v in sums.items()},all_points=points)
    for steps in (4,8,16):
        for level in (0,1):
            if (level,steps) not in field_results or (level+1,steps) not in field_results:continue
            for kind in ('native','independent'):
                pairs.append(dict(series='space',steps=steps,levels=[level,level+1],kind=kind,**delta(field_results[(level,steps)][kind],field_results[(level+1,steps)][kind])))
    for level in (0,1,2):
        for steps,next_steps in ((4,8),(8,16)):
            if (level,steps) not in field_results or (level,next_steps) not in field_results:continue
            for kind in ('native','independent'):
                pairs.append(dict(series='load',level=level,steps=[steps,next_steps],kind=kind,**delta(field_results[(level,steps)][kind],field_results[(level,next_steps)][kind])))
    write(out/'field-refinement.json',dict(passed=False,b02_complete=False,pairs=pairs))
    write(out/'index.json',dict(passed=False,b02_complete=False,cases=results,high_precision_scope='All 20 states for 768/4, original residual maximum and loaded endpoint for every other case. Original residual and field comparison cover every state.',large_6144_executed=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('output',type=Path);a=p.parse_args();analyze(a.formal.resolve(),a.output.resolve())
