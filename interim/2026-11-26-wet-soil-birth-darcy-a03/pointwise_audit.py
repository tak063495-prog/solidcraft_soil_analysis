"""Every-state difference tables for the already freshly solved pair; no solve."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from fields import CATEGORIES


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True);a=parser.parse_args()
    result=dict(observation_only_not_new_solve=True,a03_complete=False,numerical_acceptance=False,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases={})
    for name in ('hydrostatic-h1_5','hydrostatic-h2_0'):
        folder=a.run/name;sets={}
        for kind,left,right in (('accepted','native-fields.jsonl','independent-fields.jsonl'),('prepared','native-prepared.jsonl','independent-prepared.jsonl')):
            rows=[]
            with (folder/left).open() as f,(folder/right).open() as g:
                for x,y in itertools.zip_longest(f,g):
                    assert x is not None and y is not None
                    x,y=map(json.loads,(x,y));assert x['row']==y['row'] and x['fields'].keys()==y['fields'].keys()
                    stats={}
                    for key,value in x['fields'].items():
                        u,v=np.asarray(value),np.asarray(y['fields'][key])
                        assert u.shape==v.shape and np.isfinite(u).all() and np.isfinite(v).all()
                        if key in CATEGORIES:
                            stats[key]=dict(components=int(u.size),mismatches=int(np.sum(u!=v)))
                            assert stats[key]['mismatches']==0
                        else:
                            delta=u-v;maximum=float(np.max(np.abs(delta),initial=0.));absolute=float(np.linalg.norm(delta))
                            stats[key]=dict(components=int(u.size),maximum_absolute=maximum,absolute_l2=absolute,
                                            reference_l2=float(np.linalg.norm(v)),relative_l2=absolute/max(float(np.linalg.norm(v)),1e-12))
                            if key in ('u','stress','traction'):assert stats[key]['relative_l2']<1e-7
                            if key in ('internal','applied_load','effective_internal','pore_internal','contact_internal','gravity_load','surcharge_load','exterior_water_load','joint_water_load','holding_load'):
                                assert maximum<1e-10
                    rows.append(dict(row=x['row'],all_fields=stats))
            sets[kind]=rows
        result['cases'][name]=sets
    with (a.run/'pointwise-comparison.json').open('x',encoding='utf-8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    print('all accepted/prepared states checked; original gates retained; no all-field acceptance invented')
if __name__=='__main__':main()
