"""Compare the newly solved series on one identical coarsest physical-time grid."""
import argparse
import math
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from study import fields,read,write,exact,compare,sink_file


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--steps',type=int,nargs='+',default=[4,8,16]);a=p.parse_args();out=a.run.resolve();scratch=out/'common-time-grid';scratch.mkdir();report={}
    assert len(a.steps)==3 and a.steps[1]==2*a.steps[0] and a.steps[2]==2*a.steps[1]
    for label in ('VG','saturated'):
        grid=[r for r,f in fields(out/f'{label}-n{a.steps[0]}'/'native-fields.jsonl') if r['phase']=='transient'];keys={exact(r) for r in grid};assert len(keys)==4*a.steps[0]
        for kind in ('native','independent'):
            paths=[]
            for n in a.steps:
                path=scratch/f'{label}-{kind}-n{n}.jsonl';stream,sink=sink_file(path)
                try:
                    selected=[(r,f) for r,f in fields(out/f'{label}-n{n}'/(kind+'-fields.jsonl')) if exact(r) in keys]
                    assert len(selected)==len(grid)
                    for row,f in selected:sink(row,f)
                finally:stream.close()
                paths.append(path)
            pairs=[compare(paths[i],paths[i+1]) for i in range(2)];orders={}
            for k,v in pairs[0]['all_fields'].items():
                if 'difference_squared' not in v:continue
                d1=v['difference_squared']**.5;d2=pairs[1]['all_fields'][k]['difference_squared']**.5
                orders[k]=dict(coarse_to_medium_absolute_l2=d1,medium_to_fine_absolute_l2=d2,observed_order=math.log2(d1/d2) if d1>0 and d2>0 else None)
            report[label+'/'+kind]=dict(common_physical_times=grid,interpolation=False,pairs=pairs,orders=orders,
                numerical_acceptance=False,note='Fixed space; transient pair updates only. No Richardson extrapolated acceptance or JOINT time-step claim.')
    write(out/'common-grid-convergence.json',report)
if __name__=='__main__':main()
