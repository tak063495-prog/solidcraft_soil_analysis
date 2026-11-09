"""Expose auxiliary nodal holding at the force-path corner, without solving."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def audit(evidence):
    result={}
    for p in sorted(evidence.glob('*/controls.json')):
        c=json.loads(p.read_text());h=np.array(c['holding_kn']);second=np.array(c['vertical_holding_kn']);R=np.array(c['transform'])
        bounded=np.minimum(np.maximum(second,np.minimum(h,0.)),np.maximum(h,0.))
        extra=second-bounded;free=np.any(R!=0.,axis=1)
        result[p.parent.name]=dict(source_controls_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            coordinate=.5,first_group_release=1.,second_group_release=0.,
            force_decomposition=c['binding'].get('force_decomposition','unbounded_mean_force'),
            original_holding_kn=h.tolist(),corner_holding_kn=second.tolist(),
            extra_beyond_original_signed_component_bounds_kn=extra.tolist(),
            maximum_extra_magnitude_on_rows_participating_in_free_transform_kn=float(np.abs(extra[free]).max()),
            extra_generalized_force_L2_kn=float(np.linalg.norm(R.T@extra)))
    return dict(passed=False,physical_stability_certified=False,corner_measurements=result,
                scope='Auxiliary force measurement at s=.5; no material fields generated; no new acceptance threshold.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.write_text(json.dumps(audit(a.evidence),indent=2)+'\n')
