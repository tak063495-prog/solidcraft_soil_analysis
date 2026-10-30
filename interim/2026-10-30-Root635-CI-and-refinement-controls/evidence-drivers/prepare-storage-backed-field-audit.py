"""Change only result IO in the existing full seed/load/history diagnostic."""
from pathlib import Path
import ast,hashlib,json
p=Path(__file__).resolve().parent;root=p.parents[1]
out=p/'storage-backed-field-audit-v1';out.mkdir()
old=root/'runs/native-mcc-compact-pit-fixed-surface-refinement-20261029-v3/audit-fields.py'
source=old.read_text('utf-8')
anchor="analysis=p/'workflow/analysis';"
assert source.count(anchor)==1
scientific_tail=source[source.index(anchor):]
old_io="read=lambda path:json.loads(path.read_text(encoding='utf-8-sig'))\nsha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()"
new_io="""import tempfile
from export_solidcraft_view import read_result
_source_cache=tempfile.TemporaryDirectory(prefix='soil-fem-seed-load-source-')
def read(path):
 if path.name=='stage-results.json':return read_result(path,cache=Path(_source_cache.name)/'parts')
 return json.loads(path.read_text(encoding='utf-8-sig'))
def sha(path):
 with path.open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()"""
assert source.count(old_io)==1
updated=source.replace(old_io,new_io,1)
ast.parse(updated)
assert updated[updated.index(anchor):]==scientific_tail
with (out/'audit-fields.py').open('x',encoding='utf-8') as f:f.write(updated)
with (out/'original-audit-fields.py').open('x',encoding='utf-8') as f:f.write(source)
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
report=dict(prepared=True,original_driver=str(old),original_driver_sha256=sha(old),
    storage_driver_sha256=sha(out/'audit-fields.py'),
    all_original_numerical_diagnostic_statements_byte_identical=True,
    scientific_tail_sha256=hashlib.sha256(scientific_tail.encode()).hexdigest(),
    only_result_json_loading_and_whole_file_hash_io_changed=True,
    full_authoritative_result_must_still_be_hashed=True,
    no_changed_samples_thresholds_initial_stress_pc_e_removed_history_load_or_surface_checks=True,
    original_case8061_full_report_exactness_required_before34350_use=True,
    not_yet_an_actual_adapter_validation=True,p6_complete=False)
with (out/'preparation.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
print(json.dumps(dict(prepared=True,scientific_tail_byte_exact=True,not_validated_yet=True)),flush=True)
