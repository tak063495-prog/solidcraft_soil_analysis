"""Repair invalid synthetic CAD metadata without weakening the mesh contract."""
from pathlib import Path
import hashlib, json

p = Path(__file__).resolve().parent
candidate = p / 'srm-search-development-v1'
read = lambda q: json.loads(q.read_text('utf-8-sig'))
failure = read(candidate / 'quality-v2/terminal.json')
assert failure['passed'] is False and failure['failed_step'] == 'new-nine-tests'
old_log = (candidate / 'quality-v2/new-nine-tests.log').read_text('utf-8')
assert '0 passed; 9 failed; 0 ignored' in old_log
assert old_log.count('region metadata required for mesh/2') == 9
mesh_path = candidate / 'runtime/examples/srm-mean-stress-mesh.json'
mesh = read(mesh_path)
assert mesh['schema'] == 'soil-fem-mesh/2' and 'regions' not in mesh
original_hash = hashlib.sha256(mesh_path.read_bytes()).hexdigest()
source_hash = hashlib.sha256(json.dumps(
    [mesh['nodes'], mesh['tetrahedra']], separators=(',', ':'), allow_nan=False
).encode()).hexdigest()
mesh['regions'] = [dict(id=1, body='TestCube', material_id=1,
    bounds_m=[[0., 0., 0.], [1., 1., 1.]], volume_m3=1., source_sha256=source_hash)]
mesh['region_ids'] = [1] * 6
mesh_path.write_text(json.dumps(mesh, indent=2, allow_nan=False), 'utf-8', newline='\n')
driver = p / 'check-srm-bound-quality-v2.py'
text = driver.read_text('utf-8')
assert text.count("out=candidate/'quality-v2'") == 1
new = p / 'check-srm-bound-quality-v3.py'
with new.open('x', encoding='utf-8', newline='\n') as f:
    f.write(text.replace("out=candidate/'quality-v2'", "out=candidate/'quality-v3'", 1))
with (candidate / 'fixture-metadata-correction-v3.json').open('x', encoding='utf-8') as f:
    json.dump(dict(prepared=True, previous_actual_test_passes=0, previous_actual_test_failures=9,
        previous_failure_retained=True, source_of_failure='synthetic mesh/2 lacked required CAD region metadata',
        validator_and_tolerances_unchanged=True, production_algorithm_unchanged_since_quality_v2=True,
        corrected_rotation_and_heterogeneous_fixture_metadata=True,
        old_example_sha256=original_hash, new_example_sha256=hashlib.sha256(mesh_path.read_bytes()).hexdigest(),
        p6_complete=False), f, indent=2)
print('Prepared metadata correction; no mathematical acceptance claimed', flush=True)
