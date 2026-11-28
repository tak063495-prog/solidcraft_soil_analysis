"""Locate oversized native frozen partitions; replay examples are not new equilibria."""
import argparse
from collections import Counter
from pathlib import Path
from math import ceil
import numpy as np
import c01_cause_diagnostics as diag


def attribute(saved, counts):
    old,_ = diag.dependencies()
    mesh,ground,case=[diag.read(saved/name) for name in ('mesh-h0.json','ground.json','stages.json')]
    system=old.operators(mesh,ground,case,'cut').system
    result={}
    for n in (16,32,64,128,256):
        observed=diag.read(counts/f'n{n}/counts.json')
        mapping={(r['stage'],r['progress']):np.asarray(r['contact_partition_counts']) for r in observed['rows']}
        before=dict(jump=np.zeros((32,3)),traction=np.zeros((32,3)),plastic_dissipation=np.zeros(32))
        phases={s:Counter() for s in diag.STAGES};examples=[]
        for row,f in old.states(saved/f'g{n}/native'):
            actual=mapping[(row['stage'],row['progress'])]
            lower=np.array([max(1,ceil(diag.norm(v)/5e-5)) for v in f['jump']-before['jump']])
            dd=f['plastic_dissipation']-before['plastic_dissipation']
            p=phases[row['stage']]
            p['point_increments']+=32
            p['actual_above_endpoint_lower_bound']+=int(np.count_nonzero(actual>lower))
            p['actual_below_endpoint_lower_bound']+=int(np.count_nonzero(actual<lower))
            p['positive_dissipation_point_increments']+=int(np.count_nonzero(dd>0))
            p['extra_partition_and_positive_dissipation']+=int(np.count_nonzero((actual>lower)&(dd>0)))
            p['maximum_actual_count']=max(p['maximum_actual_count'],int(max(actual)))
            p['maximum_endpoint_lower_bound']=max(p['maximum_endpoint_lower_bound'],int(max(lower)))
            if len(examples)<8:
                for point in np.flatnonzero((actual>lower)&(dd>0)):
                    i=int(point);args=[before['jump'][i],before['traction'][i],f['jump'][i],system.normal[i],system.kn[i],system.kt[i],system.mu[i]]
                    a=diag.local_replay(*args,int(actual[i]))
                    b=diag.local_replay(*args,int(lower[i]))
                    examples.append(dict(stage=row['stage'],progress=row['progress'],point=i,actual_count=int(actual[i]),endpoint_lower_bound=int(lower[i]),stored_dissipation_increment=float(dd[i]),actual_partition_replay_dissipation=a[2],minimum_partition_replay_dissipation=b[2],minimum_vs_actual_traction_difference=diag.norm(np.asarray(a[0])-b[0])))
                    if len(examples)==8:break
            before=f
        result[str(n)]=dict(phases={k:dict(v) for k,v in phases.items()},first_extra_partition_plastic_examples=examples)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--saved',type=Path,required=True)
    p.add_argument('--counts',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    diag.write(a.output/'partition-attribution.json',dict(schema='c01-native-partition-attribution/1',passed=False,c01_decision='未達',new_equilibrium_solved=False,levels=attribute(a.saved,a.counts)))


if __name__=='__main__':
    main()
