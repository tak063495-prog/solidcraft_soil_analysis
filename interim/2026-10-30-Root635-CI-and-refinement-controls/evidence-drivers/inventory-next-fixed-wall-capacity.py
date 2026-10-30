"""Inventory the larger mixed system before allocating any numerical job."""
from pathlib import Path
import hashlib,json,math
p=Path(__file__).resolve().parent;root=p.parents[1]
input_root=p/'fixed-wall-finer-input-v1'
out=p/'fixed-wall-finer-capacity-preflight-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
assert read(input_root/'terminal.json')['passed']
assert read(input_root/'wall-spatial-topology-proof.json')['wall_mesh_equivalent_within_tolerance']
analysis=input_root/'workflow/analysis'
mp=analysis/'ground/mesh.json';cp=analysis/'stages.json'
mesh=read(mp);case=read(cp)
assert case['schema']=='soil-fem-projected-mcc-stages/1'
assert not case.get('initially_inactive_element_sets')
assert case['settings']['maximum_free_dofs']==100000
assert case['settings']['linear_solver']['method']=='sparse_gmres_ilu0'
assert {i for tet in mesh['tetrahedra10'] for i in tet}==set(range(len(mesh['nodes'])))
constrained=set()
for row in case['boundary']['constraints']:
    direction=row['direction']
    assert len(direction)==3 and all(x in (0.,1.) for x in direction) and sum(direction)==1.
    axis=direction.index(1.)
    assert row['value_m']==0.
    constrained.update((i,axis) for i in mesh['node_sets'][row['node_set']])
free_u=3*len(mesh['nodes'])-len(constrained)
pressure=len({(n,material) for tet,material in zip(mesh['tetrahedra'],mesh['material_ids'],strict=True) for n in tet})
mixed_primary=free_u+2*pressure
assert all(n['match_soil_node'] for n in case['beam_system']['nodes'])
candidate=root/'runs/mcc-mixed-assembly-development-20261029-v1'
quality=read(candidate/'quality-v1/source-snapshot.json')
layout=candidate/'runtime/src/stages_mixed.rs';linear=candidate/'runtime/src/linalg.rs'
assert sha(layout)==quality['src/stages_mixed.rs'] and sha(linear)==quality['src/linalg.rs']
assert 'let count = nt + dofs.count + pressure.len();' in layout.read_text('utf-8')
assert 'nt = pressure.len();' in layout.read_text('utf-8')
assert 'if self.is_dense() { 2000 } else { 300_000 }' in linear.read_text('utf-8')
assert 100000<mixed_primary<200000<300000
preview=json.loads(json.dumps(case));preview['settings']['maximum_free_dofs']=200000
restored=json.loads(json.dumps(preview));restored['settings']['maximum_free_dofs']=100000
assert restored==case
with (out/'stages-budget200k.json').open('x',encoding='utf-8') as f:json.dump(preview,f,indent=2,allow_nan=False)
report=dict(passed=True,static_inventory_only=True,nodes=len(mesh['nodes']),elements=len(mesh['tetrahedra']),
    all_mesh_nodes_initially_used=True,unique_axis_constraint_components=len(constrained),
    free_soil_displacement_dofs=free_u,material_split_P1_pressure_dofs=pressure,
    independent_projected_volume_dofs=pressure,mixed_primary_system_dimension=mixed_primary,
    installed_structural_rotations_condensed_and_not_counted_as_soil_primary=True,
    original_configured_budget=100000,original_budget_insufficient=True,
    preview_configured_budget=200000,unchanged_sparse_solver_hard_dimension_limit=300000,
    only_maximum_free_dofs_configuration_changed=True,
    all_physics_tolerances_increments_newton_and_linear_options_byte_equivalent_after_restoring_budget=True,
    original_prepared_case_unchanged=True,original_prepared_case_sha256=sha(cp),
    original_mesh_sha256=sha(mp),preview_case_sha256=sha(out/'stages-budget200k.json'),
    F2_layout_source_sha256=sha(layout),F2_linear_source_sha256=sha(linear),
    actual_allocation_resource_gate_and_solver_job_registration_and_full_solve_still_required=True,
    no_new_numerical_job_created=True,p6_complete=False)
with (out/'terminal.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2,allow_nan=False)
print(json.dumps(report),flush=True)
