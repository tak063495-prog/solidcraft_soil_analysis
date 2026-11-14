"""Inspect actual saved native Newton iterates; rejected iterates cannot prove a fold."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from fold_stability import point,low_stress
from analyze_nc_fold import material_points
from nc_fold_path import native_loads
from verify_high_order import prepare
ROOT=Path(__file__).resolve().parents[1]


def run(evidence,output):
    output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(ROOT/'interim/2026-11-09-nc-pseudo-arclength'))
    from native_arc import Native
    m=Native(evidence/'prior',Path('no-solver-called'),output/'geometry-only')
    fields=json.loads((evidence/'prior/field-records.json').read_text(encoding='utf-8'))
    matrices=json.loads((evidence/'prior/assembled-matrices.json').read_text(encoding='utf-8'))
    plane=np.asarray(json.loads((evidence/'native-coarse/samples.json').read_text())[0]['tangent'])
    final,holding=native_loads(m.mesh,m.ground,m.case,m.checkpoint['run']['continuation']['origin'],m.active)
    q,offset=prepare(m.mesh);diagnostics={};stress={}
    for prefix in ('last-accepted','terminal-rejected'):
        item=max((r for r in fields if r['label'].startswith(prefix+':')),key=lambda r:r['event']['iteration'])
        label=item['label'];event=item['event'];k=np.asarray(matrices[label+':fd:1e-09']);ad=np.asarray(matrices[label+':ad'])
        d=point(k,m.R.T@holding,plane)
        diagnostics[label]=dict(progress=event['progress'],iteration=event['iteration'],relative_force_residual=event['relative_force_residual'],
            accepted_increment=prefix=='last-accepted',use_for_fold_verdict=False,
            source='saved native map finite-difference width 1e-9 at actual Newton iterate',diagnostic=d,
            relative_matrix_difference_to_implicit_tangent=float(np.linalg.norm(k-ad)/np.linalg.norm(ad)),
            rejected_state_never_used_as_checkpoint=True)
        stress[label]=low_stress(material_points(m.mesh,m.ground,q,offset,np.asarray(m.active),event['states'],True,m.rb,np.asarray(d['right_weak_mode'])))
        stress[label].update(progress=event['progress'],iteration=event['iteration'],committed=False if prefix=='terminal-rejected' else True,
                             relative_force_residual=event['relative_force_residual'],contacts=event['joint_state'])
    (output/'native-newton-matrices.json').write_text(json.dumps(dict(diagnostic_only=True,fold_verdict='undetermined',iterates=diagnostics),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (output/'failed-low-stress.json').write_text(json.dumps(dict(diagnostic_only=True,iterates=stress),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return diagnostics,stress


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();d,s=run(a.evidence.resolve(),a.output.resolve());print(json.dumps({k:dict(p_min=s[k]['active']['pressure_kpa']['minimum'],K_min=s[k]['active']['bulk_modulus_kpa']['minimum'],sigma_min=v['diagnostic']['stiffness']['minimum']) for k,v in d.items()},indent=2))
