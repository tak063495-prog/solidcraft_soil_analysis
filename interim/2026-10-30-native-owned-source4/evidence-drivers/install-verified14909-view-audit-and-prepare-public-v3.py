"""Place the actual completed audit at its existing Native package contract path."""
from pathlib import Path
import ast, hashlib, json, shutil
p=Path(__file__).resolve().parent
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
assert not (p/'pause-requested').exists()
producer=p/'actual14909-blocked-v2'
assert read(producer/'terminal.json')['passed']
assert read(producer/'worker-passed.json')['all5_all_active_GP_VTU_comparison']
source=producer/'view-verification.json'
manifest=producer/'workflow/view/view-manifest.json'
audit, view=read(source),read(manifest)
assert audit['passed'] and audit['all_referenced_raw_points_compared_with_existing_VTU']
assert audit['raw_point_storage']=='authoritative_source'
assert audit['view_manifest_sha256']==sha(manifest)
assert audit['input_sha256']==view['input_sha256']
assert len(audit['frames'])==len(view['frames'])==5
assert [f['raw_points'] for f in audit['frames']]==[402543,402543,402543,402543,396900]
assert source.stat().st_size < 256*(1<<10)
target=producer/'workflow/view-verification.json'
assert not target.exists()
out=p/'actual14909-native-package-repair-v1'
out.mkdir()
before={str(q):sha(q) for q in [source,manifest,producer/'terminal.json',producer/'worker-passed.json']}
pending=target.with_suffix('.pending.json')
shutil.copyfile(source,pending)
assert sha(pending)==sha(source)
pending.replace(target)
assert sha(target)==sha(source) and all(sha(Path(q))==h for q,h in before.items())
old=p/'public14909-reader-v2'
assert read(old/'data-phase-terminal.json')['exit_code']==1
assert read(old/'actual-reader-terminal.json')['exit_code']!=0
assert not (old/'reader-receipt.json').exists()
failure=dict(original_session=99168,actual_compiled_public_reader_started=True,
    original_data_phase_exit_code=1,error="analysis.soil.results: file not found (os error2)",
    cause='Actual completed full VTU audit was written to producer/view-verification.json; existing Native contract reads workflow/view-verification.json',
    no_Native_check_or_mathematical_threshold_removed=True,no_result_or_audit_content_changed=True,
    exact_original_compiled_Rust_probe_reused=True,
    preserved_failed_records_sha256={name:sha(old/name) for name in ['data-phase.log','actual-reader.log','actual-reader-terminal.json','data-phase-terminal.json']})
(out/'retained-first-public-data-failure.json').write_text(json.dumps(failure,indent=2)+'\n','utf-8')
proof=dict(passed=True,only_actual_verified_audit_copied_to_existing_canonical_contract=True,
    original_report=str(source),canonical_report=str(target),report_sha256=sha(source),
    actual_all5_all_active_GP_VTU_audit_passed_before_copy=True,
    input_sha256=audit['input_sha256'],view_manifest_sha256=sha(manifest),
    original_content_bindings_unchanged=before,synthetic_or_modified_audit=False,
    no_source_or_Rust_or_solver_or_display_limit_or_gate_relaxed=True,p6_complete=False)
(out/'terminal.json').write_text(json.dumps(proof,indent=2)+'\n','utf-8')
retry=p/'public14909-reader-v3'
retry.mkdir()
for name in ['source-point-reader-probe.rs','compiled-probe-binding.json','public-source-point-reader.exe']:
    shutil.copyfile(old/name,retry/name)
    assert sha(old/name)==sha(retry/name)
original_script=p/'read-real14909-v2.py'
text=original_script.read_text('utf-8')
old_literal="out=p/'public14909-reader-v2'"
new_literal="out=p/'public14909-reader-v3'"
assert text.count(old_literal)==1
updated=text.replace(old_literal,new_literal,1)
assert ast.dump(ast.parse(updated.replace(new_literal,old_literal,1)))==ast.dump(ast.parse(text))
script=p/'read-real14909-v3.py'
script.write_text(updated,'utf-8',newline='\n')
(retry/'preparation.json').write_text(json.dumps(dict(prepared=True,
    original_data_comparison_driver_AST_identical_except_output_path=True,
    source_code_and_compiled_Rust_probe_and_engine_and_serde_byte_exact=True,
    all525_original_raw_node_derived_plate_and_inactive_and_display_settings_checks_retained=True,
    first_failed_data_phase_and_packaging_repair_recorded=True,
    script_sha256=sha(script),prior_script_sha256=sha(original_script),p6_complete=False),indent=2),'utf-8')
print(json.dumps(proof),flush=True)
