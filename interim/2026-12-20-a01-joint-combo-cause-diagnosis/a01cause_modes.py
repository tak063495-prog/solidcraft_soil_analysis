"""Mode localization and nominal stiffness ratios on captured systems only."""
import argparse
import json
from pathlib import Path
import numpy as np
from a01cause_binding import fixtures,original_reference
from a01cause_linear import native_events,matrix,assess,rigid_fit
from a01cause_driver import write
from a01cause_metrics import norm,dot


def main():
    p=argparse.ArgumentParser();p.add_argument('--formal',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    folder=a.formal/'stacked';mesh,ground,case=(fixtures.read(folder/(key+'.json')) for key in ('mesh','ground','case'))
    model=original_reference.Model(mesh,ground,case);event=native_events(folder/'native/trace.jsonl')[-1]
    original=event['original_tangent'];transform=np.zeros((len(original['dof_rows']),len(event['rhs'])))
    for i,row in enumerate(original['dof_rows']):
        for j,c in row:transform[i,j]=c
    skeleton=transform.T@model.k@transform;s=np.linalg.svd(skeleton,compute_uv=False)
    u,sv,vh=np.linalg.svd(matrix(event));projection=dot(u[:,-1],event['rhs'])
    young=[x['mohr_coulomb']['young_kpa'] for x in ground['materials']]
    length=float(np.ptp(np.asarray(mesh['nodes'])[:,2]))/3
    records=[json.loads(line) for line in (folder/'independent/trace.jsonl').read_text(encoding='utf-8').splitlines()]
    birth=[r for r in records if r['kind']=='linear' and r['stage']=='birth']
    independent=[]
    for e in (birth[0],birth[-1]):
        _,_,right=np.linalg.svd(e['matrix']);independent.append(dict(iteration=e['iteration'],
            rigid_fits=rigid_fit(mesh,np.asarray(e['soil_t']),right[-3:]),assessment=assess(e['matrix'],e['rhs'])))
    write(a.output,dict(passed=False,observational_only=True,acceptance_tolerance_changed=False,
        continuum_only_reduced_dimension=len(s),continuum_only_numerical_rank=int(np.sum(s>len(s)*np.finfo(float).eps*s[0])),
        continuum_only_smallest_singular_values=s[-8:].tolist(),
        nominal_block_height_m=length,young_kpa=young,kn_kpa_per_m=model.system.kn.tolist(),kt_kpa_per_m=model.system.kt.tolist(),
        nominal_kn_L_over_E_min=float(model.system.kn.min()*length/max(young)),nominal_kn_L_over_E_max=float(model.system.kn.max()*length/min(young)),
        actual_smallest_mode_left_load_projection=projection,actual_smallest_mode_required_coefficient=projection/sv[-1],
        actual_smallest_mode_rigid_fits=rigid_fit(mesh,transform,[vh[-1]]),independent_birth_mode_localization=independent,
        physical_instability_proved=False))

if __name__=='__main__':main()
