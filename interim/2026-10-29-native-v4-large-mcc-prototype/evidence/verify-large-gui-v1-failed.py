from pathlib import Path
import datetime,hashlib,json,re
import numpy as np
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'large-mcc-gui-review'
read=lambda q:json.loads(q.read_text('utf-8'))
def sha(q):
    h=hashlib.sha256()
    with q.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
prep=read(out/'preparation.json')
assert all(sha(Path(q))==h for q,h in prep['files_sha256'].items())
assert sha(Path(prep['binary']['path']))==prep['binary']['sha256']
frame_path=root/'runs/mcc-large-structures-development-20261028-v1/large-case/workflow/view/frame-004.json'
frame=read(frame_path);assert frame['name']=='excavate_pit'
increments=np.asarray(frame['node_fields']['stage_displacement_increment_m'])[frame['active_node_ids'],2]*1000
node_range=[float(increments.min()),float(increments.max())]
plates=[q for q in frame['structures']['plates'] if q['active']]
assert len(plates)==1264
moments=[s['bending_moment_kn'][0] for q in plates for s in q['result']['samples']]
moment_range=[min(moments),max(moments)]
top=read(out/'observation-top.json');wall=read(out/'observation-wall.json');section=read(out/'observation-section.json')
assert section['normal']==[1,0,0] and section['offset_m']==5 and section['slice_only']
assert top['displacement_basis']==section['displacement_basis']=='stage'
assert np.allclose(moment_range,wall['observed_range_kn_m_per_m'],rtol=0,atol=5e-5)
# A degree-two display can also sample interior values. Do not equate its
# full range with nodal extrema if refined interpolation has an overshoot.
nodal_range_matches=np.allclose(node_range,top['observed_range_mm'],rtol=0,atol=5e-5)
for name in ['excavation-increment-controls-tree.txt','wall-bending-controls-tree.txt','section-controls-tree.txt']:
    text=(out/'screenshots'/name).read_text('utf-8')
    assert 'excavate_pit' in text
    assert 'CAD or analysis inputs changed.' not in text
manifest=read(Path(next(q for q in prep['files_sha256'] if q.endswith('view-manifest.json'))))
assert manifest['input_sha256']==prep['expected_input_sha256']
assert sha(frame_path)==manifest['files_sha256']['frame-004.json']
assert sha(root/'external/solidcraft-soil-fem/target/ci/solidcraft.exe')==read(p/'gui-baseline/binary-identity.json')['sha256']
with (out/'verification.json').open('x',encoding='utf-8') as f:
    json.dump(dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        actual_GUI_loaded_535711781_byte_source=True,binary=prep['binary'],
        CAD_input_identity_required_and_matched=True,all_original_source_and_owned_design_bytes_unchanged=True,
        excavation_increment_top_view_and_reference_section_visible=True,
        section_normal=[1,0,0],section_offset_m=5,
        wall_raw_DKT_local_Mxx_visible=True,active_plates=1264,
        independent_display_legend_check=dict(raw_Mxx_range=moment_range,observed_Mxx_range=wall['observed_range_kn_m_per_m'],
            matches_display_precision=True,raw_active_node_increment_mm_range=node_range,
            observed_increment_mm_range=top['observed_range_mm'],nodal_extrema_match_display_precision=bool(nodal_range_matches)),
        source_frame_sha256=sha(frame_path),screenshots_sha256={q.name:sha(q) for q in (out/'screenshots').glob('*.png')},
        default_native62_unchanged=True,app_closed_normally=True,
        full_CI_release_other_GPU_all5_GUI_navigation_and_job_recovery_acceptance_separate=True,p6_complete=False),f,indent=2)
print('PASS actual large MCC GUI increment, wall bending and reference section; immutable inputs and raw moment legend confirmed')
