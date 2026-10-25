"""Diagnostic field comparison at an accepted initialization fraction, not full-run acceptance."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def compare(mesh_file,native_file,coarse_file,fine_file,output):
    mesh,native,coarse,fine=map(read,(mesh_file,native_file,coarse_file,fine_file))
    geo=native['geostatic'];actual=geo['last_committed'];fraction=geo['progress']
    assert fraction==fine['load_fraction']==coarse['load_fraction']
    assert actual['active_elements']==fine['active_elements'] and not native['complete']
    for reference in (coarse,fine):
        for name,digest in reference['input_and_source_sha256'].items():assert sha(Path(name))==digest,name
    state=np.asarray(fine['material_state']).reshape(-1,14);seed=np.asarray(fine['seed_material_state']).reshape(-1,14)
    start=geo['seed']['material_states'];assert len(state)==len(actual['material_states'])==len(start)
    assert np.allclose(seed[:,:6],[s['effective_stress_compression_positive_kpa'] for s in start],rtol=0,atol=1e-12)
    fields={}
    def metric(a,b):
        a=np.asarray(a);b=np.asarray(b);assert a.shape==b.shape
        return dict(relative_l2=float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30)),maximum_absolute=float(abs(a-b).max()))
    for name,expected,observed in (
        ('stress_kpa',state[:,:6],[s['effective_stress_compression_positive_kpa'] for s in actual['material_states']]),
        ('pc_kpa',np.exp(state[:,6]),[s['preconsolidation_pressure_kpa'] for s in actual['material_states']]),
        ('void_ratio',np.exp(state[:,7])-1,[s['void_ratio'] for s in actual['material_states']]),
        ('plastic_strain',state[:,8:],[s['plastic_strain_compression_positive'] for s in actual['material_states']])):
        fields[name]=metric(observed,expected)
    p=state[:,:3].mean(axis=1);pn=np.array([s['effective_stress_compression_positive_kpa'] for s in actual['material_states']])[:,:3].mean(axis=1)
    plastic=np.array([s['plastic_strain_compression_positive'] for s in actual['material_states']])
    inventory=[];cursor=0;seed_dof={(r['node_index'],r['material_id']):r['compression_positive_kpa'] for r in geo['seed']['mixed_mcc']['pressure_dofs']}
    seed_projection=[]
    for e,tet in enumerate(mesh['tetrahedra']):
        q=ref.volume(mesh,e);nq=len(q['weight']);l=np.column_stack((1-q['natural'].sum(axis=1),q['natural']))
        for j in range(nq):
            pt=native['integration_points'][cursor+j];assert pt['element_index']==e and np.allclose(pt['position_m'],q['position'][j],rtol=0,atol=1e-8)
        seed_projection.extend(l@np.array([seed_dof[(node,mesh['material_ids'][e])] for node in tet]))
        inventory.extend(range(cursor,cursor+nq));cursor+=nq
    assert cursor==len(state)
    point=int(np.argmax(np.linalg.norm(plastic-state[:,8:],axis=1)))
    report=dict(schema='soil-fem-mcc-initial-fraction-comparison/1',diagnostic_only=True,goal_complete=False,p6_complete=False,
        full_initialization_accepted=False,load_fraction=fraction,native_run_complete=native['complete'],
        common_seed_stresses_verified=True,integration_layout_verified=True,
        independent_reference_steps=[coarse['reference_steps'],fine['reference_steps']],
        reference_displacement_refinement=metric(coarse['displacement_m'],fine['displacement_m']),
        native_displacement=metric(actual['displacement_m'],fine['displacement_m']),fields=fields,
        reference_maximum_strain=fine['maximum_engineering_strain_from_seed'],
        native_maximum_material_strain=geo['maximum_absolute_engineering_strain_from_seed'],
        native_plastic_points=int(np.count_nonzero(np.linalg.norm(plastic,axis=1)>1e-10)),
        independent_plastic_points=int(np.count_nonzero(np.linalg.norm(state[:,8:],axis=1)>1e-10)),
        initial_pressure_projection_maximum_error_kpa=float(abs(seed[:,:3].mean(axis=1)-seed_projection).max()),
        largest_plastic_difference=dict(integration_point=point,position_m=native['integration_points'][point]['position_m'],
            reference_mean_kpa=float(p[point]),native_mean_kpa=float(pn[point]),
            reference_pc_kpa=float(np.exp(state[point,6])),native_state=actual['material_states'][point]),
        input_and_source_sha256={str(p):sha(p) for p in (mesh_file,native_file,coarse_file,fine_file,Path(__file__),Path(ref.__file__))})
    with output.open('x',encoding='utf-8') as f:json.dump(report,f,indent=2,allow_nan=False)
    print(json.dumps({k:report[k] for k in ('load_fraction','reference_displacement_refinement','native_displacement','native_plastic_points','independent_plastic_points','initial_pressure_projection_maximum_error_kpa')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('native',type=Path);p.add_argument('coarse',type=Path);p.add_argument('fine',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();compare(*(getattr(a,k).resolve() for k in ('mesh','native','coarse','fine','output')))
