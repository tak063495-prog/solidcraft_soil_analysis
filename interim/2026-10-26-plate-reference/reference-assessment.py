"""Measure unchanged elastic rates near the actual input K0 state."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

folder=Path(__file__).resolve().parent;root=folder.parents[1]
sys.path.insert(0,str(root/'tools'))
import freeze_sources
import mcc_rate_fem as rate
import mcc_projected_sparse_reference as fresh
import mcc_projected_factor_reference as cached
import mcc_projected_matrixfree_reference as matrixfree

def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

freeze_sources.freeze(folder);sources=freeze_sources.manifest()
inputs=root/'runs/projected-k0-reset-reference-20261026-v1/inputs'
mesh,ground,case,identity=[read(inputs/n) for n in ('mesh.json','ground.json','case.json','mesh-identity.json')]
active=np.ones(len(mesh['tetrahedra']),bool)
state=rate.seed(mesh,ground,identity).reshape(-1,14);plastic=np.zeros(len(state),bool)
target=case['stages'][0]['target_pressure_kpa']
force=rate.physical(mesh,ground,active,target,case.get('pressure_on_original_surface_only',False))[0].reshape(-1)-rate.physical(mesh,ground,active,ground['initialization']['surface_surcharge_kpa'],case.get('pressure_on_original_surface_only',False))[0].reshape(-1)
states=[]
for i in range(12):
    candidate=state.copy();candidate[:,:3]*=1+i*1e-5
    assert rate.terms(candidate,np.array([[ground['materials'][int(m)-1]['mcc']['critical_state_m']**2,
        ground['materials'][int(m)-1]['mcc']['lambda'],ground['materials'][int(m)-1]['mcc']['kappa'],
        ground['materials'][int(m)-1]['mcc']['poisson']] for m in np.repeat(mesh['material_ids'],27)]))[0].max()<0.
    states.append(candidate)
from scipy.sparse.linalg import splu
class AlternateFresh(fresh.SparseProjectedSystem):
    def factorize(self,matrix):return splu(matrix,permc_spec='MMD_AT_PLUS_A')
systems=[fresh.SparseProjectedSystem(mesh,ground,active,case),cached.ReusedFactorProjectedSystem(mesh,ground,active,case),matrixfree.MatrixFreeFactorProjectedSystem(mesh,ground,active,case),AlternateFresh(mesh,ground,active,case)]
outputs=[];times=[]
for system in systems:
    system.rate(state,force,plastic) # Warm both paths before timing.
    start=time.perf_counter();values=[system.rate(candidate,force,plastic) for candidate in states]
    times.append(time.perf_counter()-start);outputs.append(values)
comparisons=[]
for index,name in [(1,'reused_LU'),(2,'matrix_free_reused_LU'),(3,'alternate_fresh_LU_ordering')]:
    maximum_relative=0.;strict=True;max_component=0.;group_errors={}
    for ordinary,other in zip(outputs[0],outputs[index],strict=True):
        for a,b in zip(ordinary,other,strict=True):
            strict &= bool(np.allclose(a,b,rtol=1e-8,atol=1e-10))
            maximum_relative=max(maximum_relative,float(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-12)))
            max_component=max(max_component,float(np.max(np.abs(a-b))))
        for key,a,b in [('displacement_rate',ordinary[0],other[0]),('stress_rate',ordinary[1][:,:6],other[1][:,:6]),
            ('log_volume_rate',ordinary[1][:,7],other[1][:,7])]:
            group_errors[key]=max(group_errors.get(key,0.),float(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-30)))
    comparisons.append(dict(name=name,strict_componentwise_gate_satisfied=strict,
        maximum_relative_array_difference=maximum_relative,maximum_component_difference=max_component,relative_physical_group_errors=group_errors))
scaling=[]
for name,values in zip(('fresh','reused','matrix_free','alternate_fresh'),outputs,strict=True):
    errors=[]
    for i,value in enumerate(values):
        alpha=1+i*1e-5;start=values[0]
        errors.append(dict(displacement=float(np.linalg.norm(value[0]-start[0]/alpha)/max(np.linalg.norm(start[0]/alpha),1e-30)),
            stress=float(np.linalg.norm(value[1][:,:6]-start[1][:,:6])/max(np.linalg.norm(start[1][:,:6]),1e-30)),
            strict_componentwise_stress_scaling=bool(np.allclose(value[1][:,:6],start[1][:,:6],rtol=1e-8,atol=1e-10))))
    scaling.append(dict(name=name,maximum_relative_displacement_scaling_error=max(v['displacement'] for v in errors),
        maximum_relative_stress_scaling_error=max(v['stress'] for v in errors),strict_componentwise_stress_scaling_all_satisfied=all(v['strict_componentwise_stress_scaling'] for v in errors)))
for system in systems:
    assert system.maximum_true_linear_residual<=1e-10
    assert max(system.maximum_block_residuals.values())<=1e-10
assert freeze_sources.manifest()==sources
report=dict(schema='soil-fem-independent-reused-factor-rate-assessment/1',original_equation_residual_gates_passed=True,
    scope='12 elastic-rate perturbations near the analytical K0 state of the actual CAD input; NOT a completed nonlinear construction reference',
    nodes=len(mesh['nodes']),tet10=len(mesh['tetrahedra']),material_points=len(state),
    fresh_seconds=times[0],reused_seconds=times[1],measured_warm_rate_speedup=times[0]/times[1],
    matrix_free_seconds=times[2],measured_warm_matrix_free_speedup=times[0]/times[2],
    comparisons=comparisons,analytical_uniform_stiffness_scaling=scaling,
    fresh=systems[0].solver_diagnostics(),reused=systems[1].solver_diagnostics(),
    matrix_free=systems[2].solver_diagnostics(),
    input_sha256={n:sha(inputs/n) for n in ('mesh.json','ground.json','case.json','mesh-identity.json')},
    source_manifest_sha256=sha(folder/'source-snapshot.json'),source_files=len(sources),
    strict_componentwise_rate_equivalence_accepted=comparisons[0]['strict_componentwise_gate_satisfied'],
    all_nonlinear_fields_accepted=False,goal_complete=False,p6_complete=False)
(folder/'assessment.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps(report),flush=True)
