"""Measure a user-PC global sequence, without fitting an acceptance threshold."""
import argparse
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('pc_step_measure',HERE/'measure_steps.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def summarize(folder,steps,h,maximum,output):
    assert steps==sorted(set(steps)) and len(steps)>=3 and all(b==2*a for a,b in zip(steps,steps[1:]))
    pairs={};orders={};entities={}
    for kind in ('native','independent'):
        def path(g):return folder/f'h{h}-g{g}-m{maximum:g}-cut'/kind
        for a,b in zip(steps,steps[1:]):
            key=f'{kind}/{a}__{b}';pairs[key]=m.old.compare(path(a),path(b));entities[key]=m.entity_differences(path(a),path(b))
        for gs in zip(steps,steps[1:],steps[2:]):
            row=m.old.richardson([path(g) for g in gs])
            for value in row['fields'].values():
                for k in list(value):
                    if k.startswith('richardson_correction'):del value[k]
            orders[kind+'/'+ '__'.join(map(str,gs))]=row
    output.mkdir(parents=True,exist_ok=False)
    m.write(output/'measurements.json',dict(passed=False,steps=steps,pairs=pairs,observed_orders=orders,all_entity_differences=entities,
             continuous_solution_certified=False,reason='PC measurements only. No event alignment, extrapolated field or new acceptance tolerance. C-02 pointwise spatial mapping is outside this global-sequence summary.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--steps',type=int,nargs='+',required=True);p.add_argument('--h',type=int,default=0);p.add_argument('--maximum',type=float,default=5e-5)
    a=p.parse_args();summarize(a.input.resolve(),a.steps,a.h,a.maximum,a.output.resolve())
