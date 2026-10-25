"""Independently audit the material/kinematic small-strain domain at initialization.

Incomplete runs are allowed for diagnosis. Completing this audit never turns
a numerical failure into a physical failure or an accepted equilibrium state.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def audit(mesh_file,result_file,output):
    mesh=json.loads(mesh_file.read_text(encoding='utf-8-sig'))
    result=json.loads(result_file.read_text(encoding='utf-8-sig'))
    geo=result['geostatic'];state=geo['last_committed'];mixed=state['mixed_mcc']
    active=np.asarray(state['active_elements']);u=np.asarray(state['displacement_m'])
    theta=np.asarray(mixed['cumulative_volume_strain_compression_positive']);layout=result['integration_points']
    material=[];kinematic=[];point_ids=[];cursor=0
    for e in range(len(mesh['tetrahedra'])):
        q=ref.volume(mesh,e);nq=len(q['weight'])
        if active[e]:
            kin=-np.einsum('qai,i->qa',q['b'],u[q['nodes']].reshape(-1))
            mat=kin.copy();mat[:,:3]+=(theta[cursor:cursor+nq]-kin[:,:3].sum(axis=1))[:,None]/3
            for local in range(nq):
                assert layout[cursor+local]['element_index']==e
                assert np.allclose(layout[cursor+local]['position_m'],q['position'][local],rtol=0.,atol=1e-8)
            material.extend(mat);kinematic.extend(kin);point_ids.extend(range(cursor,cursor+nq))
        cursor+=nq
    assert cursor==len(layout)==len(theta)
    material=np.asarray(material);kinematic=np.asarray(kinematic)
    assert np.isfinite(material).all() and np.isfinite(kinematic).all()
    row,component=np.unravel_index(np.argmax(abs(material)),material.shape);point=point_ids[row]
    maximum=float(abs(material).max());kin_max=float(abs(kinematic).max())
    limit=result['ground']['initialization']['geostatic']['maximum_initial_engineering_strain']
    report=dict(schema='soil-fem-mixed-initial-strain-audit/1',audit_completed=True,diagnostic_only=True,
        goal_complete=False,p6_complete=False,physical_failure_inferred=False,
        native_run_complete=result['complete'],geostatic_complete=geo['complete'],geostatic_progress=geo['progress'],
        configured_maximum_initial_engineering_strain=limit,maximum_material_engineering_strain_from_seed=maximum,
        maximum_kinematic_engineering_strain_from_seed=kin_max,
        cached_geostatic_maximum_strain=geo['maximum_absolute_engineering_strain_from_seed'],
        initial_state_within_configured_strain_domain=bool(max(maximum,kin_max)<=limit),
        largest_material_strain=dict(integration_point=point,component=int(component),position_m=layout[point]['position_m'],
            material_engineering_strain=material[row].tolist(),kinematic_engineering_strain=kinematic[row].tolist(),
            cumulative_material_volume_strain=float(theta[point]),state=state['material_states'][point]),
        input_and_source_sha256={str(p):sha(p) for p in (mesh_file,result_file,Path(__file__),Path(ref.__file__))})
    with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2,allow_nan=False)
    print(json.dumps({k:report[k] for k in ('native_run_complete','geostatic_progress','configured_maximum_initial_engineering_strain',
        'maximum_material_engineering_strain_from_seed','maximum_kinematic_engineering_strain_from_seed','initial_state_within_configured_strain_domain')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('result',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();audit(a.mesh.resolve(),a.result.resolve(),a.output.resolve())
