"""C-01 style exact saved-gap interval labels applied to C-02 common points."""
import argparse
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import c02diag_driver as d
from c02diag_mapping import ContactPlan
from c02diag_labels import catalog,masks,budget,direction_budget

FIELDS=('plastic_slip','traction','plastic_dissipation','numerical_dissipation')


def run(source,out):
    out.mkdir(parents=True,exist_ok=False);c=d.legacy()
    meshes,ground,case,systems,positions,roots,weights=d.common(source,c)
    plans={(h,m):ContactPlan(meshes[h],systems[h],h,positions,roots,m) for h in range(3) for m in d.METHODS}
    routes=[(h,g) for h in (0,1) for g in (32,64,128,256)]+[(2,g) for g in (32,64)]
    histories={};catalogs={};mapped={}
    for h,g in routes:
        rows,fields=d.load_contact(c,source/f'h{h}-g{g}'/'native');histories[h,g]=(rows,fields)
        gap=plans[h,'nearest'].history(fields['gap']);catalogs[h,g]=catalog(rows,gap)
        for m in d.METHODS:mapped[h,g,m]={k:plans[h,m].history(fields[k]) for k in FIELDS}
    result=dict(passed=False,new_global_equilibrium_solved=False,labels_are_telescoping_accounting_not_causal_effects=True,local_replay_used_as_equilibrium=False,source='saved GP histories from sealed 2026-11-19',spatial={},time_spatial_change={})
    saved={}
    for g in (32,64):
        mask_pairs=[masks([catalogs[h,g],catalogs[h+1,g]],c.STAGES,g) for h in (0,1)]
        for m in d.METHODS:
            components={}
            for h,(mask,duration) in enumerate(mask_pairs):
                key=f'fixed{g}-{m}-h{h}{h+1}'
                row=dict(label_point_times=int(mask.sum()),largest_open_set_duration_difference=float(np.max(duration)),fields={})
                for k in FIELDS:
                    delta=mapped[h,g,m][k]-mapped[h+1,g,m][k]
                    summary,part=budget(delta,mask,c.STAGES,g);row['fields'][k]=summary
                    components[h,k]=part
                    if k in ('plastic_slip','traction'):
                        for name,array in part.items():saved[key+'-'+k+'-'+name]=array
                result['spatial'][key]=row
            result['spatial'][f'fixed{g}-{m}-direction_cross_terms']={k:direction_budget(components[0,k],components[1,k]) for k in FIELDS}
    for (ha,hb),ga,gb in [((0,1),32,64),((1,2),32,64),((0,1),64,256),((0,1),128,256)]:
        route=[(ha,ga),(hb,ga),(ha,gb),(hb,gb)]
        mask,duration=masks([catalogs[x] for x in route],c.STAGES,ga)
        for m in d.METHODS:
            values=[d.sampled(histories[h,g][0],mapped[h,g,m],ga)[1] for h,g in route]
            key=f'h{ha}{hb}-g{ga}-{gb}-{m}'
            result['time_spatial_change'][key]=dict(label_point_times=int(mask.sum()),largest_open_set_duration_difference=float(np.max(duration)),fields={})
            for k in FIELDS:
                delta=(values[2][k]-values[3][k])-(values[0][k]-values[1][k])
                summary,part=budget(delta,mask,c.STAGES,ga)
                result['time_spatial_change'][key]['fields'][k]=summary
                if m=='cubic' and k=='plastic_slip':
                    for name,array in part.items():saved[key+'-'+k+'-'+name]=array
    d.write(out/'label-budget.json',result)
    d.write(out/'open-intervals.json',{f'h{h}-g{g}':catalogs[h,g] for h,g in routes})
    np.savez_compressed(out/'label-components.npz',**saved)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',required=True,type=Path);parser.add_argument('--output',required=True,type=Path);a=parser.parse_args();run(a.source.resolve(),a.output.resolve())
