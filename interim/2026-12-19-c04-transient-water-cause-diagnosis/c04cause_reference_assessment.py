"""Level differences on exactly the same input coordinates; no extrapolated pass."""
import argparse
import math
from pathlib import Path
import c04cause_metrics as metrics
from c04cause_observe import common_pair
from c04cause_driver import write
FIELDS=('u','traction','plastic_slip','accumulated','plastic_dissipation','numerical_dissipation','endpoint_work')


def assess(series):
    common=set((r['stage'],r['progress']) for r in series[0][0])
    for rows,_ in series[1:]:common&={(r['stage'],r['progress']) for r in rows}
    ordered=[(r['stage'],r['progress']) for r in series[0][0] if (r['stage'],r['progress']) in common]
    arrays=[]
    for rows,fields in series:
        index={(r['stage'],r['progress']):i for i,r in enumerate(rows)}
        arrays.append({k:v[[index[p] for p in ordered]] for k,v in fields.items()})
    report={}
    for key in FIELDS:
        differences=[arrays[i][key]-arrays[i+1][key] for i in range(len(arrays)-1)]
        norms=[metrics.norm(v) for v in differences]
        report[key]=dict(adjacent_absolute_l2=norms,
          adjacent_relative_l2=[metrics.metric(arrays[i][key],arrays[i+1][key])['relative_l2'] for i in range(len(arrays)-1)],
          estimated_orders=[math.log2(a/b) if a>0 and b>0 else None for a,b in zip(norms,norms[1:])],
          adjacent_difference_cosines=[metrics.cosine(a,b) for a,b in zip(differences,differences[1:])])
    return dict(common_accepted_coordinates=len(ordered),interpolation=False,fields=report,
                accepted_continuum_reference=False,richardson_correction_applied=False,new_acceptance_gate=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--controls',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();base=a.baseline;out=a.controls;report={}
    for kind in ('native','independent'):
        report[kind+'-water-8-16-32']=assess([metrics.load(base/f'{kind}-w{n}/fields.jsonl') for n in (8,16,32)])
    for kind in ('native','adaptive'):
        report[kind+'-mechanical-32-64-128']=assess([metrics.load(out/(f'{kind}-w32' if n==32 else f'{kind}-m{n}-w32')/'fields.jsonl') for n in (32,64,128)])
    report['independent-fixed-contact-1-2-4-8']=assess([metrics.load(out/f'fixed{n}-w32/fields.jsonl') for n in (1,2,4,8)])
    report['passed']=False;report['finest_continuum_reliability_established']=False
    write(a.report,report)

if __name__=='__main__':main()
