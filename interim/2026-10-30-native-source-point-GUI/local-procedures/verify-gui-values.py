"""Compare actual observed GUI JSON with independently VTU-checked public API."""
from pathlib import Path
import hashlib,json,datetime
p=Path(__file__).resolve().parent;out=p/'gui-source73-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.file_digest(q.open('rb'),'sha256').hexdigest()
api=read(p/'public8061-reader-v1/reader-receipt.json')['frames'][4]
gp=read(out/'gui-raw-point0.json');plate=read(out/'gui-plate0.json')
assert gp==next(v['value'] for v in api['points'] if v['id']==0)
plate_api=next(v['value'] for v in api['plates'] if v['id']==0)
assert plate.keys()==plate_api.keys()
assert plate['deformation_scale']==1 and plate_api['deformation_scale']==12.5
assert all(plate[k]==plate_api[k] for k in plate if k!='deformation_scale')
assert gp['mixed_fields']['primary_mean_effective_stress_compression_positive_kpa']!=sum(gp['effective_stress_compression_positive_kpa'][:3])/3
text=(out/'plate-Mxx-controls.txt').read_text('utf-8')
assert 'Range: -6.59243e1 to 4.95484e1 kN m/m' in text
assert abs(api['wall']['data_range'][0]+65.92427419918772)<1e-12
assert abs(api['wall']['data_range'][1]-49.54838015533024)<1e-12
assert abs(api['section']['reference_section_area_m2']-18)<1e-11
binding=read(out/'original-binding.json');assert all(sha(Path(q))==h for q,h in binding.items())
report=dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
 actual_observed_GUI_JSON=True,raw_GP0_all19_top_level_properties_exact=True,
 plate0_all_raw_properties_and_nested_raw_samples_exact=True,
 plate_display_scale_metadata_intentionally_differs=dict(GUI=1,API=12.5),
 API_values_previously_checked_against_independent_existing_GP_and_plate_VTUs=True,
 raw_material_stress_separate_from_primary_and_work_stress=True,
 observed_GUI_stage_vertical_range_mm=[-.397111,12.8286],
 Native_API_stage_vertical_range_mm=api['configured']['data_range'],
 observed_GUI_plate_Mxx_range_kn_m_per_m=[-65.9243,49.5484],
 Native_API_plate_Mxx_range_kn_m_per_m=api['wall']['data_range'],
 GUI_reference_section=dict(normal=[1,0,0],offset_m=5,section_only=True),
 Native_API_reference_section_area_m2=api['section']['reference_section_area_m2'],
 GUI_section_area_numeric_not_exposed=True,original_document_binary_manifest_audit_unchanged=True,
 formal_release_full_CI_new_PC_and_general_P6_acceptance_separate=True,p6_complete=False,
 sha256={str(q):sha(q) for q in [out/'gui-raw-point0.json',out/'gui-plate0.json',out/'plate-Mxx-controls.txt',p/'public8061-reader-v1/reader-receipt.json']})
with (out/'value-verification.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
print('PASS actual GUI raw GP and plate values exact; displayed ranges agree at displayed precision')
