"""Summarize the three completed full-field cases at the original criteria."""
from pathlib import Path
import hashlib,json,math
p=Path(__file__).resolve().parent;study=p/'plane-strain-slope-comparison-v1'
out=study/'completed-controls-comparison-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def norm(rows,weights):return math.sqrt(math.fsum(w*x*x for row,w in zip(rows,weights,strict=True) for x in row))
def delta(a,b):
    assert len(a)==len(b) and all(len(x)==len(y) for x,y in zip(a,b))
    return [[x-y for x,y in zip(u,v)] for u,v in zip(a,b)]
names=['xy1-ny1','xy1-ny2','xy1-ny1-global64'];data={};assessment={};binding={}
protocol=read(study/'protocol.json');binding[str(study/'protocol.json')]=sha(study/'protocol.json')
for name in names:
    folder=study/name;field=folder/'field-comparison-v1'
    assert read(folder/'solve-v1/terminal.json')['passed'] and read(field/'full-independent-audit-terminal.json')['passed']
    assessment[name]=read(field/'terminal.json');data[name]=read(field/'field-samples.json')
    assert assessment[name]['original_predeclared_limits']==protocol['predeclared_acceptance']
    assert assessment[name]['protocol_sha256']==sha(study/'protocol.json')
    assert assessment[name]['gates']['actual_rejected_history_unchanged']
    for q in [field/'terminal.json',field/'field-samples.json',field/'full-independent-audit-terminal.json',
        folder/'mesh.json',folder/'ground.json',folder/'stages.json']:binding[str(q)]=sha(q)
assert read(study/'xy1-ny1/stages.json')==read(study/'xy1-ny2/stages.json')
assert read(study/'xy1-ny1/ground.json')==read(study/'xy1-ny2/ground.json')
for name in ['xy1-ny1','xy1-ny2']:
    mesh=read(study/name/'mesh.json')
    assert max(n[1] for n in mesh['nodes'])-min(n[1] for n in mesh['nodes'])==1.
base=data['xy1-ny1'];ref=base['source2d_u_xz_m'];stress=base['source2d_prism_mean_stress_xyz_xy_yz_xz_kpa']
weights=base['original200_quad_prism_volumes_m3'];depths=['0.0','0.5','1.0']
assert len(ref)==661 and len(stress)==len(weights)==200
unorm=norm(ref,[1.]*661);umax=max(math.hypot(*r) for r in ref);snorm=norm(stress,weights)
recomputed={};pairs={}
for name,s in data.items():
    assert s['source2d_u_xz_m']==ref and s['source2d_prism_mean_stress_xyz_xy_yz_xz_kpa']==stress
    assert set(s['native3d_u_xz_m'])==set(depths)
    assert all(len(s['native3d_u_xz_m'][d])==661 for d in depths)
    assert all(abs(a-b)<1e-10 for a,b in zip(weights,s['original200_quad_prism_volumes_m3']))
    error={d:norm(delta(s['native3d_u_xz_m'][d],ref),[1.]*661)/unorm for d in depths}
    assert all(math.isclose(error[d],assessment[name]['metrics']['displacement_sample_relative_L2_per_depth'][d],rel_tol=1e-12,abs_tol=1e-14) for d in depths)
    recomputed[name]=error
    if name!='xy1-ny1':
        pairs[name]=dict(displacement_change_over_fixed_2d_reference_L2_per_depth={d:norm(delta(s['native3d_u_xz_m'][d],base['native3d_u_xz_m'][d]),[1.]*661)/unorm for d in depths},
            maximum_displacement_change_over_fixed_2d_reference_max=max(math.hypot(*v)/umax for d in depths for v in delta(s['native3d_u_xz_m'][d],base['native3d_u_xz_m'][d])),
            prism_mean_stress_change_over_fixed_2d_reference_L2=norm(delta(s['native3d_prism_mean_stress_xyz_xy_yz_xz_kpa'],base['native3d_prism_mean_stress_xyz_xy_yz_xz_kpa']),weights)/snorm)
assert all(sha(Path(q))==h for q,h in binding.items())
report=dict(passed=True,diagnostic_comparison_completed=True,all_full_mechanical_history_raw_VTU_audits_passed=True,
    all661_nodes_each_of3_depths_and_all200_prisms_compared=True,
    original_reference_acceptance_by_case={n:assessment[n]['passed'] for n in names},
    independently_recomputed_2d_displacement_relative_L2=recomputed,pairwise_changes_from_baseline=pairs,
    thickness_width_fixed_at_one_metre_not_widened=True,
    reported_depth_variation={n:assessment[n]['metrics']['depth_displacement_variation_over_max'] for n in names},
    original5percent_and10percent_limits_unchanged=True,
    section_refinement_and_material_increment_control_pending=True,
    no_asymptotic_convergence_or_cause_or_physical_Fs_upper_claim=True,p6_complete=False)
for name,value in [('terminal.json',report),('source-binding.json',binding)]:
    with (out/name).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False)
print(json.dumps(report),flush=True)
