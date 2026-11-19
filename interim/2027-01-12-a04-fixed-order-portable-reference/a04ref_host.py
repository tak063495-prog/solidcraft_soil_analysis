"""Fresh Linux host baselines from fixture inputs and exact-coordinate field differences."""
import argparse
import json
from pathlib import Path
from a04ref_matrix import driver
import numpy as np


def variants():
    for saturated in (False,True):
        label = 'saturated' if saturated else 'VG'
        for n in (4,8,16):
            yield 'soil/'+label+'-n'+str(n),n,3,saturated,False
        yield label+'-reverse-keys',4,3,saturated,True
        for order in (5,7):
            yield label+'-quadrature'+str(order),8,order,saturated,False


def run(out):
    if np.finfo(np.longdouble).nmant<=52:
        raise ValueError('this comparison requires an extended host baseline')
    out.mkdir(parents=True,exist_ok=False)
    study=driver.bind()['study']
    for label,n,order,saturated,reverse in variants():
        folder=out/label
        folder.mkdir(parents=True)
        for k,v in zip(('mesh','ground','case'),study.fixture(n,saturated)):
            driver.write(folder/(k+'.json'),driver.reverse_keys(v) if reverse else v)
        driver.independent(folder,order)
        print(label,flush=True)


def compare(host,portable):
    reports={}
    for label,n,order,saturated,reverse in variants():
        a,b=host/label,portable/label
        for name in ('mesh','ground','case'):
            assert driver.read(a/(name+'.json'))==driver.read(b/(name+'.json'))
        maxima={}
        with (a/'independent-fields.jsonl').open() as stream_a,(b/'independent-fields.jsonl').open() as stream_b:
            for count,(x,y) in enumerate(zip(map(json.loads,stream_a),map(json.loads,stream_b),strict=True),1):
                assert x['row']==y['row'] and x['fields'].keys()==y['fields'].keys()
                for k in x['fields']:
                    va,vb=np.array(x['fields'][k]),np.array(y['fields'][k])
                    assert va.shape==vb.shape
                    diff=np.abs(va.astype(float)-vb.astype(float))
                    maxima[k]=max(maxima.get(k,0.),float(diff.max(initial=0.)))
        trace=driver.read(b/'a04ref-line-search.json')['iterations']
        reports[label]=dict(states=count,maximum_absolute=maxima,
            accepted_hydraulic_steps=sum(not row['trials'] for row in trace),
            newton_updates=sum(bool(row['trials']) for row in trace),
            line_search_trials=sum(len(row['trials']) for row in trace),
            maximum_halvings=max((len(row['trials'])-1 for row in trace),default=0),
            maximum_final_merit=max(float(row['merit']) for row in trace if not row['trials']),
            maximum_true_linear_residual=max(row.get('linear',0.) for row in trace))
    return reports


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--host',required=True,type=Path)
    p.add_argument('--compare',type=Path)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.compare:
        driver.write(a.output,compare(a.host,a.compare))
    else:
        run(a.host)
