"""Recompute explicit field/refinement gates; never infer full P6 completion."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from compare_mcc_rate_fem import compare


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def assess(run,mesh_path,lifecycle,output,candidate=None,candidate_case=None):
    if output.exists():raise ValueError('refusing to overwrite a previous assessment')
    mesh=read(mesh_path);references=[read(run/f'slope-cycle-reference-{n}.json') for n in (25,50,100)]
    candidate=candidate or run/'slope-cycle-rust-finer.json'
    candidate_case=candidate_case or run/'slope-cycle-finer-case.json'
    natives=[read(run/'slope-cycle-rust.json'),read(candidate)]
    comparisons=[compare(references[-1],native,mesh) for native in natives]
    gates=[]
    def gate(name,condition,observed=None,limit=None):
        gates.append(dict(name=name,passed=bool(condition),observed=observed,limit=limit))
    gate('geometry',len(mesh['nodes'])==459 and len(mesh['tetrahedra10'])==224)
    coarse_case=read(run/'slope-cycle-case.json');fine_case=read(candidate_case)
    coarse_physics=json.loads(json.dumps(coarse_case));fine_physics=json.loads(json.dumps(fine_case))
    for key in ('material_maximum_substep','initial_increment','minimum_increment'):
        coarse_physics['settings'].pop(key);fine_physics['settings'].pop(key)
    gate('same_physics_after_numerical_refinement',coarse_physics==fine_physics)
    gate('candidate_numerical_resolution',0<fine_case['settings']['material_maximum_substep']<=1e-6
        and 0<fine_case['settings']['initial_increment']<=1/32
        and 0<fine_case['settings']['minimum_increment']<=fine_case['settings']['initial_increment'])
    for native in natives:
        gate('completed_four_frames',native['complete'] and len(native['frames'])==4)
        gate('same_mesh_identity',native['provenance']['mesh_sha256']==read(run/'mesh-identity.json')['mesh_sha256'])
    lifecycle_report=read(lifecycle/'verification.json')
    gate('cycle_rollback_and_restarts',lifecycle_report['passed'] and lifecycle_report['plastic_reload_rejection']
         and len(lifecycle_report['separate_process_restarts'])==5)
    for name,digest in lifecycle_report['evidence_sha256'].items():
        gate('lifecycle_evidence:'+name,sha(lifecycle/name)==digest)
    for i,row in enumerate(comparisons[-1]['frames']):
        prefix=row['name']+':'
        gate(prefix+'plastic_point_count',row['independent_plastic_points']==row['rust_plastic_points']==(7 if i==0 else 188))
        for key,limit in (('displacement_relative_all_nodes',3e-5),('displacement_maximum_component_error_m',1e-6),
                          ('display_displacement_maximum_component_error_m',1e-6),('independent_relative_force_residual',1e-9)):
            gate(prefix+key,row[key]<limit,row[key],limit)
        if i:
            key='stage_displacement_maximum_component_error_m';gate(prefix+key,row[key]<1e-6,row[key],1e-6)
        for field,key,limit in (('stress','relative_all_points',1e-5),('pc','maximum_absolute_all_points',1e-4),
                                ('void_ratio','maximum_absolute_all_points',5e-6),('plastic_strain','maximum_absolute_all_points',1e-6)):
            value=row['fields'][field][key];gate(prefix+field+':'+key,value<limit,value,limit)
        gate(prefix+'rust_refinement_reduces_displacement_error',
             row['displacement_relative_all_nodes']<comparisons[0]['frames'][i]['displacement_relative_all_nodes'])
    refinement=[]
    for i in range(4):
        arrays=[np.asarray(data['frames'][i]['displacement_m']) for data in references]
        changes=[float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30)) for a,b in zip(arrays,arrays[1:])]
        gate(f'reference_frame_{i}_refinement',changes[1]<changes[0] and changes[1]<1e-5,changes,1e-5)
        refinement.append(dict(name=references[-1]['frames'][i]['name'],relative_displacement_changes=changes))
    removed=mesh['excavation_sets']['Pit']
    for reference in references:
        before=np.asarray(reference['frames'][-2]['material_state']);after=np.asarray(reference['frames'][-1]['material_state'])
        gate('removed_reference_histories_frozen',np.array_equal(before[removed],after[removed]))
    for native in natives:
        before=native['frames'][-2]['state']['material_states'];after=native['frames'][-1]['state']['material_states']
        gate('removed_native_histories_frozen',all(before[e*27+q]==after[e*27+q] for e in removed for q in range(27)))
    files=[mesh_path,lifecycle/'verification.json',Path(__file__),Path(__file__).with_name('compare_mcc_rate_fem.py')]
    files += [run/f'slope-cycle-reference-{n}.json' for n in (25,50,100)]
    files += [run/name for name in ('slope-cycle-rust.json','slope-cycle-case.json','mesh-identity.json')]+[candidate,candidate_case]
    result=dict(schema='soil-fem-mcc-cycle-refinement-assessment/1',passed=all(g['passed'] for g in gates),
        goal_complete=False,p6_complete=False,scope='Dry small-strain 5 kPa Tet10 slope load/unload/local cut; numerical/field gates only',
        gates=gates,reference_refinement=refinement,comparisons=comparisons,
        sha256={str(p):sha(p) for p in files},
        limitations='Requires separate source/input provenance and archived evidence before scoped acceptance; no spatial convergence or general P6 claim')
    with output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
    print(json.dumps(dict(passed=result['passed'],failed_gates=[g for g in gates if not g['passed']],reference_refinement=refinement)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('mesh',type=Path)
    p.add_argument('lifecycle',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--candidate',type=Path);p.add_argument('--candidate-case',type=Path);a=p.parse_args()
    assess(*(v.resolve() for v in (a.run,a.mesh,a.lifecycle,a.output)),
           a.candidate.resolve() if a.candidate else None,a.candidate_case.resolve() if a.candidate_case else None)
