"""Check all base-tracked bytes, existing seals/records and the production CLI."""
import hashlib
import json
from pathlib import Path
import subprocess
out = Path(__file__).resolve().parent
root = out.parents[1]
base = 'b020b492ea091288dbfffac0514c570898f8c408'
raw = subprocess.check_output(['git', 'ls-tree', '-r', '-z', base], cwd=root)
rows = []
for entry in raw.split(b'\0'):
    if not entry: continue
    meta, raw_name = entry.split(b'\t', 1); name = raw_name.decode(); blob = meta.decode().split()[2]
    path = root/name
    data = path.read_bytes()
    actual = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    rows.append(dict(path=name, base_blob=blob, current_blob=actual, unchanged=blob == actual))
changed = [r['path'] for r in rows if not r['unchanged']]
assert set(changed) == {'.github/workflows/fem-ci.yml', 'tools/verify_wet_structures.py'}
prior = [r for r in rows if r['path'].startswith('interim/')]
assert all(r['unchanged'] for r in prior)
for prefix in ('src/', 'examples/', 'docs/'):
    assert all(r['unchanged'] for r in rows if r['path'].startswith(prefix))
prior_binary = json.loads((root/'interim/2026-11-08-nc-full-pivot-and-fixed-prescribed/binary-and-source-sha256.json').read_text())['production']['target/release/soil-fem-3d']
binary = hashlib.sha256((root/'target/release/soil-fem-3d').read_bytes()).hexdigest()
assert binary == prior_binary
report = dict(base_commit=base, base_tracked_files=len(rows), unchanged_base_files=sum(r['unchanged'] for r in rows),
              changed_existing_files=changed, old_interim_files=len(prior), old_interim_all_byte_identical=True,
              existing_docs_examples_numerical_rust_all_byte_identical=True,
              production_binary_sha256=binary, production_binary_matches_base_evidence=True,
              existing_seals_reassigned=False, linux_seals_reassigned_to_windows=False,
              rows=rows)
(out/'base-preservation.json').write_text(json.dumps(report, indent=2)+'\n')
print({k: v for k, v in report.items() if k != 'rows'})
