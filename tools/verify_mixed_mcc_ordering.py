"""Compare physical states at the same accepted mixed increment across ordering.

This is a solver-ordering check, not an independent material/field acceptance.
Elapsed time is reported without a hardware-dependent pass criterion.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))


def verify(before,after,output):
    a,b=read(before)['run'],read(after)['run'];x,y=a['last_committed'],b['last_committed']
    ca,cb=copy.deepcopy(a['case']),copy.deepcopy(b['case'])
    oa,ob=ca['settings'].pop('linear_solver'),cb['settings'].pop('linear_solver')
    assert ca==cb and a['ground']==b['ground'],'ordering comparison changed physical or material-integration inputs'
    assert a['provenance']['mesh_sha256']==b['provenance']['mesh_sha256'],'different meshes'
    assert a['continuation']==b['continuation'],'different construction state'
    assert a['integration_points']==b['integration_points'],'different integration layout'
    assert x['active_elements']==y['active_elements']
    for key in ('applied_nodal_load_kn','excavation_release_load_kn','stage_index','stage_progress','boundary_factor','pressure_kpa'):
        if key in x or key in y:assert x[key]==y[key],key
    rows=[]
    def maximum(name,left,right,tolerance):
        left=np.asarray(left);right=np.asarray(right);assert left.shape==right.shape
        error=float(np.max(abs(left-right)));assert np.isfinite(error) and error<tolerance,name
        rows.append(dict(name=name,maximum_absolute_difference=error,limit=tolerance,passed=True))
    maximum('displacement_m',x['displacement_m'],y['displacement_m'],1e-10)
    for key,tolerance in (('effective_stress_compression_positive_kpa',1e-7),('preconsolidation_pressure_kpa',1e-7),
                          ('void_ratio',1e-10),('plastic_strain_compression_positive',1e-10)):
        maximum(key,[s[key] for s in x['material_states']],[s[key] for s in y['material_states']],tolerance)
    mx,my=x['mixed_mcc'],y['mixed_mcc']
    assert mx['tangent_method']==my['tangent_method']=='implicit_branch_fixed_partition'
    assert mx['last_increment_partitions']==my['last_increment_partitions']
    for key in ('cumulative_volume_strain_compression_positive','reference_void_ratio'):
        maximum(key,mx[key],my[key],1e-10)
    assert [(r['node_index'],r['material_id']) for r in mx['pressure_dofs']]==[(r['node_index'],r['material_id']) for r in my['pressure_dofs']]
    maximum('primary_pressure_kpa',[r['compression_positive_kpa'] for r in mx['pressure_dofs']],
            [r['compression_positive_kpa'] for r in my['pressure_dofs']],1e-7)
    numerical=[]
    for tag,run,options in (('before',a,oa),('after',b,ob)):
        solves=[s for step in run['attempts'] for s in step['linear_solves']]
        assert solves and all(step['accepted'] for step in run['attempts'])
        assert max(s['true_relative_residual'] for s in solves)<=run['case']['settings']['relative_linear_tolerance']
        numerical.append(dict(name=tag,options=options,solves=len(solves),
            methods=sorted({s['method'] for s in solves}),maximum_iterations=max(s['iterations'] for s in solves),
            total_iterations=sum(s['iterations'] for s in solves),elapsed_seconds=sum(s['elapsed_seconds'] for s in solves),
            maximum_original_relative_residual=max(s['true_relative_residual'] for s in solves)))
    assert numerical[1]['total_iterations']<numerical[0]['total_iterations']
    assert numerical[1]['methods']==['csr_gmres_ilu0_volume_first_rcm_up']
    report=dict(schema='soil-fem-mixed-mcc-ordering-comparison/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Same accepted dry mixed-MCC column increment, numerical ordering and fixed-partition history only',
        field_differences=rows,numerical=numerical,material_partitions_exact=True,
        solver_before_sha256=a['provenance']['solver_binary_sha256'],solver_after_sha256=b['provenance']['solver_binary_sha256'],
        input_and_source_sha256={str(p):sha(p) for p in (before,after,Path(__file__))},
        limitations='Different explicitly recorded Krylov widths; no global timing guarantee or independent constitutive acceptance')
    with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(dict(passed=True,field_differences=rows,numerical=numerical)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('before',type=Path);p.add_argument('after',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();verify(a.before.resolve(),a.after.resolve(),a.output.resolve())
