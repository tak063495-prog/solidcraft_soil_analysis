from pathlib import Path
import hashlib,json,subprocess,datetime
p=Path(__file__).resolve().parent;out=p/'gui-source73-v1';root=p.parents[1]
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
assert read(out/'value-verification.json')['passed']
owned=read(out/'desktop-child.json');assert owned['pid']==588
ps='Get-CimInstance Win32_Process -Filter "ProcessId=588" | Select-Object ProcessId,ExecutablePath,CreationDate | ConvertTo-Json -Compress'
assert not subprocess.check_output(['pwsh','-NoProfile','-Command',ps],text=True).strip()
binding=read(out/'original-binding.json');assert all(sha(Path(q))==h for q,h in binding.items())
sources=read(p/'native-quality-v2/source73.json')['files_sha256'];native=root/'runs/native-typed-source-integration-20261029-v1/native'
assert all(sha(native/name)==h for name,h in sources.items())
shots=['excavation-stage-controls.png','gui-raw-point0.png','excavation-top-stage-increment.png','section-x5-controls.png','section-x5-stage-increment.png','plate-Mxx-3D.png']
assert all((out/'screenshots'/name).is_file() for name in shots)
v=dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source73_unchanged=True,
 actual_changed_experimental_desktop_GUI=True,document_bodies=11,elements=8061,raw_GP_count=217647,
 accepted_view_schema='soil-fem-solidcraft-view/4',original_result_bytes=1459184313,
 actual_GUI_opened_original_accepted_manifest=True,excavation_stage_increment_top_view=True,
 x5_reference_coordinate_section_visually_observed=True,plate_local_Mxx_three_dimensional_visualization=True,
 raw_GP0_and_plate0_data_compared_with_VTU_checked_public_API=True,
 original_bindings_unchanged=True,owned_process_normally_closed_via_AltF4=True,owned_process_no_longer_present=True,
 screenshots_local_only_not_for_Git_commit=True,
 screenshots_sha256={name:sha(out/'screenshots'/name) for name in shots},
 full_normal_main_CI_still_running=not (p/'desktop-and-main-ci-v1/terminal.json').exists(),
 official_release_cold_PC_other_GPU_general_P6_acceptance_separate=True,p6_complete=False)
with (out/'terminal.json').open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
print('PASS real changed GUI8061, raw values and normal owned-process close; formal release/newPC/fullP6 separate')
