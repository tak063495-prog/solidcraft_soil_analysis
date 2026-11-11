"""Reproducible, bounded independent slope diagnostics (NumPy/SciPy only)."""
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path
import time
import numpy as np
import scipy
import lem
import srm
from compare import MODEL, UNITS, CONVENTIONS, REDUCTION, case_fingerprint

HERE=Path(__file__).resolve().parent


def dump(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n','utf-8')


def write_fields(directory,result):
    eq=result['equilibrium_lower']
    if eq is None: return None
    trial=next(t for t in result['trials'] if t['factor']==eq and t['status']=='equilibrated')
    xy=np.array(result['mesh']['xy']); tri=np.array(result['mesh']['triangles']); u=np.array(trial['displacement'])
    from compare import NODE, PLASTIC
    for filename,columns,rows in [
        ('nodes.csv',NODE,[[i+1,*p,*v] for i,(p,v) in enumerate(zip(xy,u))]),
        ('plastic.csv',PLASTIC,[[i+1,*p,*v] for i,(p,v) in enumerate(zip(xy[tri].mean(axis=1),trial['plastic_strain']))])]:
        with (directory/filename).open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f); w.writerow(columns); w.writerows(rows)
    return dict(factor=eq,gravity_fraction=1.,reference='stress_free_gravity_turn_on',
                plastic_measure='engineering_tensor_components',nodes_csv='nodes.csv',plastic_csv='plastic.csv')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--grids',type=int,nargs='+',default=[3,5,9]);p.add_argument('--slices',type=int,nargs='+',default=[16,32,64])
    p.add_argument('--h-levels',type=int,nargs='+',default=[4,6,8]);p.add_argument('--trial-seconds',type=float,default=2.)
    p.add_argument('--bounds',type=float,nargs=6,metavar=('A0','A1','B0','B1','S0','S1'))
    p.add_argument('--max-trials',type=int,default=8);p.add_argument('--gravity-steps',type=int,default=8)
    p.add_argument('--width',type=float,default=.025);p.add_argument('--lem-only',action='store_true');p.add_argument('--srm-only',action='store_true')
    args=p.parse_args()
    if min(args.grids)<2 or min(args.slices)<4 or min(args.h_levels)<1 or args.trial_seconds<=0:
        p.error('positive budgets and valid grid dimensions required')
    if args.lem_only and args.srm_only: p.error('mutually exclusive modes')
    args.output.mkdir(parents=True,exist_ok=False)
    start=time.monotonic(); case_hash=case_fingerprint()
    manifest=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                  case_sha256=case_hash,accepted=None,fs_tolerance='未定義',P6_02='未達',
                  controls=vars(args)|{'output':str(args.output)},source_sha256={q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in HERE.glob('*.py')})
    dump(args.output/'environment.json',manifest)
    if not args.srm_only:
        values=[]
        for nc in [False,True]:
            for grid in args.grids:
                for count in args.slices:
                    bounds=tuple(zip(args.bounds[::2],args.bounds[1::2])) if args.bounds else ((2.,11.5),(20.,32.),(1.,7.))
                    r=lem.search(grid,count,nc,bounds); values.append(r)
                    dump(args.output/'lem.json',values)
                    print('LEM',nc,grid,count,{k:None if v is None else v['fs'] for k,v in r['minima'].items()},flush=True)
    if not args.lem_only:
        series=[]
        for n in args.h_levels:
            r=srm.bisect_srm(n,width=args.width,max_trials=args.max_trials,seconds=args.trial_seconds,steps=args.gravity_steps)
            folder=args.output/f'srm-n{n}';folder.mkdir();dump(folder/'raw.json',r)
            fields=write_fields(folder,r)
            trials=[dict(factor=t['factor'],status=t['status'],reason=t['reason'],relative_residual=t['relative_residual'],
                         residual_tolerance=1e-7,gravity_fraction=1. if t['status']=='equilibrated' else t['completed_steps']/args.gravity_steps,
                         failure_certified=False) for t in r['trials']]
            interval=r['computational_transition']
            interchange=dict(schema_version=1,model_id=MODEL,case_sha256=case_hash,units=UNITS,conventions=CONVENTIONS,
                             reduction=REDUCTION,accepted=None,fs_tolerance=None,trials=trials,trial_count=len(trials),
                             factor_interval=interval,endpoint_statuses=[next(t['status'] for t in trials if t['factor']==f) for f in interval] if interval else None,
                             interval_kind='computational_unresolved' if interval else None,field_state=fields)
            dump(folder/'interchange.json',interchange)
            brief={k:v for k,v in r.items() if k not in ['mesh','trials']};brief['trial_states']=[[t['factor'],t['status']] for t in r['trials']]
            series.append(brief);dump(args.output/'srm-series.json',series)
            print('SRM',n,r['equilibrium_lower'],r['unresolved_upper'],r['certified_fs_interval'],flush=True)
    dump(args.output/'legacy-diagnosis.json',dict(
        original_intervals=[[1.3250000000000002,1.3531250000000001],[1.3320312500000002,1.3337890625000002]],
        original_interval_model='phi=0,15 kPa uniform stress column; analytical Fs=4/3, NOT this slope',
        against_original_analytical_4_over_3=['含まれる','含まれる'],
        against_slope_reference='判定不能: different physical problem, numerical membership would be misleading',
        original_central_L2_percent=[15.582,15.886,15.593,15.672],original_threshold=.05,
        original_discrete_accepted=[False]*4,original_continuum_accepted=None,
        reason='No matched in-plane h series or original 661-node fields solved here',accepted=None,P6_02='未達'))
    manifest['elapsed_seconds']=time.monotonic()-start;dump(args.output/'environment.json',manifest)

if __name__=='__main__':main()
