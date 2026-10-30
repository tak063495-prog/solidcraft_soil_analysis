"""Quantify the completed loading control, with the original field failure retained."""
from pathlib import Path
import hashlib,json,math
p=Path(__file__).resolve().parent
study=p/'plane-strain-slope-comparison-v1'
out=study/'completed-loading-control-comparison-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def norm(rows,weights=None):
    if weights is None:weights=[1.]*len(rows)
    assert len(rows)==len(weights)
    return math.sqrt(math.fsum(w*x*x for row,w in zip(rows,weights) for x in row))
def diff(a,b):
    assert len(a)==len(b) and all(len(x)==len(y) for x,y in zip(a,b))
    return [[x-y for x,y in zip(row,other)] for row,other in zip(a,b)]
names=['xy1-ny1','xy1-ny1-global64']
protocol=read(study/'protocol.json');bindings={str(study/'protocol.json'):sha(study/'protocol.json')}
assessments=[];samples=[]
for name in names:
    folder=study/name
    assert read(folder/'solve-v1/terminal.json')['passed']
    assert read(folder/'field-comparison-v1/full-independent-audit-terminal.json')['passed']
    a=read(folder/'field-comparison-v1/terminal.json');s=read(folder/'field-comparison-v1/field-samples.json')
    assert a['passed'] is False and a['original_predeclared_limits']==protocol['predeclared_acceptance']
    assert a['fixed_factor']==protocol['reference_fixed_factor']
    assert a['gates']['actual_rejected_history_unchanged'] and a['gates']['full_independent_audit']
    assert not a['gates']['all_661_displacements_each_depth'] and not a['gates']['maximum_displacement_error']
    assessments.append(a);samples.append(s)
    for item in ['mesh.json','ground.json','stages.json','reference-map.json',
        'solve-v1/terminal.json','field-comparison-v1/terminal.json','field-comparison-v1/field-samples.json',
        'field-comparison-v1/full-independent-audit-terminal.json','field-comparison-v1/binding.json']:
        q=folder/item;bindings[str(q)]=sha(q)
base,control=[study/n for n in names]
for name in ['mesh.json','ground.json','reference-map.json']:assert sha(base/name)==sha(control/name)
original=read(base/'stages.json');modified=read(control/'stages.json')
changes={k:modified['settings'][k] for k in original['settings'] if original['settings'][k]!=modified['settings'][k]}
assert changes==dict(initial_increment=1/64,minimum_increment=1/8192)
restored=json.loads(json.dumps(modified));restored['settings']=original['settings'];assert restored==original
a,b=samples;ref=a['source2d_u_xz_m'];stress_ref=a['source2d_prism_mean_stress_xyz_xy_yz_xz_kpa']
weights=a['original200_quad_prism_volumes_m3'];depths=['0.0','0.5','1.0']
assert ref==b['source2d_u_xz_m'] and stress_ref==b['source2d_prism_mean_stress_xyz_xy_yz_xz_kpa']
assert len(ref)==661 and len(stress_ref)==len(weights)==200
assert weights==b['original200_quad_prism_volumes_m3']
assert all(len(a['native3d_u_xz_m'][d])==len(b['native3d_u_xz_m'][d])==661 for d in depths)
unorm=norm(ref);umax=max(math.hypot(*v) for v in ref);snorm=norm(stress_ref,weights)
assert min(unorm,umax,snorm)>0
pair_u={d:norm(diff(a['native3d_u_xz_m'][d],b['native3d_u_xz_m'][d]))/unorm for d in depths}
pair_max=max(math.hypot(*v)/umax for d in depths for v in diff(a['native3d_u_xz_m'][d],b['native3d_u_xz_m'][d]))
pair_stress=norm(diff(a['native3d_prism_mean_stress_xyz_xy_yz_xz_kpa'],b['native3d_prism_mean_stress_xyz_xy_yz_xz_kpa']),weights)/snorm
assert all(math.isfinite(v) for v in list(pair_u.values())+[pair_max,pair_stress])
assert all(sha(Path(q))==h for q,h in bindings.items())
report=dict(passed=True,comparison_completed=True,only_loading_controls_changed=changes,
    same_all661_nodes_all200_prisms_material_load_boundary_mesh_and_FSS=True,
    both_full_independent_force_history_geostatic_and_raw_VTU_audits_passed=True,
    displacement_change_over_fixed_2d_reference_L2_per_depth=pair_u,
    maximum_displacement_change_over_fixed_2d_reference_max=pair_max,
    prism_mean_stress_change_over_fixed_2d_reference_L2=pair_stress,
    reference_field_acceptance_passed=False,
    original_2d_relative_L2_at_middle=[q['metrics']['displacement_sample_relative_L2_per_depth']['0.5'] for q in assessments],
    actual_cutbacks_with_history_restored=[q['metrics']['actual_geostatic_rejections'] for q in assessments],
    original_displacement_relative_L2_limit=.05,original_maximum_displacement_limit=.1,
    same_or_smaller_material_increment_not_evaluated_by_this_control=True,
    mesh_convergence_or_cause_or_physical_Fs_upper_not_inferred=True,
    pairwise_change_is_diagnostic_not_a_new_accuracy_acceptance=True,p6_complete=False)
for name,value in [('terminal.json',report),('source-binding.json',bindings)]:
    with (out/name).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False)
print(json.dumps(report),flush=True)
