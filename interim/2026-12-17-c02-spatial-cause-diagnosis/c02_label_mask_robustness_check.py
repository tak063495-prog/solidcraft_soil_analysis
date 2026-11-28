"""Check the broader four-path opening label against within-mesh time labels.

Only reads already sealed catalogs/differences. No equilibrium or local replay.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from c02diag_labels import masks,budget
STAGES=('geostatic','shear','open','recontact','reverse','reload','return')


def run(labels,out):
    out.mkdir(parents=True,exist_ok=False)
    intervals=json.loads((labels/'open-intervals.json').read_text(encoding='utf-8'))
    result=dict(passed=False,scope='plastic slip time dependence of spatial differences, cubic reconstruction',opening_labels_are_accounting_not_causal=True,comparisons={})
    with np.load(labels/'label-components.npz',allow_pickle=False) as arrays:
        for (ha,hb),ga,gb in [((0,1),32,64),((1,2),32,64),((0,1),64,256),((0,1),128,256)]:
            key=f'h{ha}{hb}-g{ga}-{gb}-cubic'
            left,dla=masks([intervals[f'h{ha}-g{ga}'],intervals[f'h{ha}-g{gb}']],STAGES,ga)
            right,dlb=masks([intervals[f'h{hb}-g{ga}'],intervals[f'h{hb}-g{gb}']],STAGES,ga)
            within=left|right
            four,_=masks([intervals[f'h{h}-g{g}'] for h,g in [(ha,ga),(hb,ga),(ha,gb),(hb,gb)]],STAGES,ga)
            total=arrays[key+'-plastic_slip-opening_label']+arrays[key+'-plastic_slip-other_path']
            strict,_=budget(total,within,STAGES,ga);broad,_=budget(total,four,STAGES,ga)
            result['comparisons'][key]=dict(within_mesh_time_labels=strict,four_path_labels=broad,within_mesh_label_point_times=int(within.sum()),four_path_label_point_times=int(four.sum()),largest_within_mesh_open_set_duration_difference=float(max(np.max(dla),np.max(dlb))))
    result['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [labels/'open-intervals.json',labels/'label-components.npz']}
    with (out/'label-mask-robustness.json').open('x',encoding='utf-8',newline='\n') as f:json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--labels',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();run(a.labels.resolve(),a.output.resolve())
