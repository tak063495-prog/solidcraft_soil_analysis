"""Compare input-only contact-discretization variants at the same physical points."""
import argparse
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import c02diag_driver as d
from c02diag_discretization import alternate_system
from c02diag_mapping import ContactPlan
from c02diag_transfer import norm,dot


def run(source,variants,out):
    out.mkdir(parents=True,exist_ok=False);c=d.legacy()
    meshes,ground,case,systems,positions,roots,weights=d.common(source,c)
    baseline={};variant_fields={};metadata={};result=dict(passed=False,new_contact_algorithm_introduced=False,input_only_global_equilibrium=True,production_default_changed=False,per_mesh={},spatial={},time_interaction={})
    configurations=[('baseline','jacobi4','tri6'),('legendre4-tri6','legendre4','tri6'),('jacobi4-corner3','jacobi4','corner3')]
    for name,rule,interpolation in configurations:
        for h in (0,1):
            system=systems[h] if name=='baseline' else alternate_system(systems[h],rule,interpolation)
            plans={m:ContactPlan(meshes[h],system,h,positions,roots,m) for m in d.METHODS}
            for g in (32,64):
                folder=source/f'h{h}-g{g}'/'independent' if name=='baseline' else variants/name/f'h{h}-g{g}'/'independent'
                rows,fields=d.load_contact(c,folder)
                if len(rows)!=7*g:raise ValueError('incomplete variant equilibrium path')
                for m in d.METHODS:variant_fields[name,h,g,m]={k:plans[m].history(v) for k,v in fields.items()}
                if name!='baseline':
                    path=variants/name/f'h{h}-g{g}';audit=d.read(path/'result.json')
                    if audit['reference_uses_native_results'] or audit['maximum_relative_force_residual']>=1e-10:raise ValueError('variant equilibrium criteria failed')
                    metadata[name,h,g]=dict(contact_gp=len(system.points),accepted_steps=len(audit['path']),maximum_relative_force_residual=audit['maximum_relative_force_residual'],rule=rule,interpolation=interpolation,layout_sha256=d.sha(path/'layout.json'),manifest_sha256=d.sha(folder/'manifest.json'))
    for name,rule,interpolation in configurations[1:]:
        for h in (0,1):
            for g in (32,64):
                key=f'{name}-h{h}-g{g}'
                result['per_mesh'][key]=dict(audit=metadata[name,h,g],methods={m:{k:d.difference(variant_fields[name,h,g,m][k],variant_fields['baseline',h,g,m][k]) for k in d.FIELDS} for m in d.METHODS})
        for g in (32,64):
            for m in d.METHODS:
                key=f'{name}-g{g}-{m}';result['spatial'][key]={}
                for k in d.FIELDS:
                    a=variant_fields[name,0,g,m][k];b=variant_fields[name,1,g,m][k]
                    x=variant_fields['baseline',0,g,m][k];y=variant_fields['baseline',1,g,m][k]
                    change=(a-b)-(x-y);scale=norm(x-y)
                    result['spatial'][key][k]=dict(baseline=d.difference(x,y),variant=d.difference(a,b),spatial_difference_change_absolute_l2=norm(change),change_over_baseline_spatial_difference=norm(change)/max(scale,1e-12),difference_cosine=dot((a-b).ravel(),(x-y).ravel())/max(norm(a-b)*scale,1e-24))
        for m in d.METHODS:
            key=f'{name}-{m}';result['time_interaction'][key]={}
            for k in d.FIELDS:
                def change(n):
                    coarse=variant_fields[n,0,32,m][k]-variant_fields[n,1,32,m][k]
                    fine=variant_fields[n,0,64,m][k][1::2]-variant_fields[n,1,64,m][k][1::2]
                    return fine-coarse
                base=change('baseline');modified=change(name)
                result['time_interaction'][key][k]=dict(baseline_time_change_absolute_l2=norm(base),variant_time_change_absolute_l2=norm(modified),difference_of_time_changes_l2=norm(modified-base))
    d.write(out/'discretization-summary.json',result)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True,type=Path);p.add_argument('--variants',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();run(a.source.resolve(),a.variants.resolve(),a.output.resolve())
