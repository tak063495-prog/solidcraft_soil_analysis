"""Verify every new archive byte, then replay sealed audit code in isolation."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((HERE/'archive-sha256.json').read_text(encoding='utf-8'))
    archive = out/'joined.zip'; digest = hashlib.sha256(); size = 0
    with archive.open('xb') as destination:
        for part in manifest['parts']:
            raw = (HERE/part['name']).read_bytes()
            if len(raw) != part['size'] or sha(raw) != part['sha256']:
                raise ValueError('changed archive part')
            destination.write(raw); digest.update(raw); size += len(raw)
    if size != manifest['archive_size'] or digest.hexdigest() != manifest['archive_sha256']:
        raise ValueError('changed joined archive')
    evidence = out/'evidence'; evidence.mkdir()
    seen = set()
    with zipfile.ZipFile(archive) as z:
        for entry in z.infolist():
            path = Path(entry.filename)
            if path.is_absolute() or '..' in path.parts or entry.is_dir() or entry.filename in seen:
                raise ValueError('unsafe/duplicate archive member')
            expected = manifest['entries'].get(entry.filename)
            if expected is None:
                raise ValueError('unlisted archive member')
            raw = z.read(entry)
            if len(raw) != expected['size'] or sha(raw) != expected['sha256']:
                raise ValueError('changed archive entry')
            target = evidence/path; target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(raw)
            seen.add(entry.filename)
    if seen != set(manifest['entries']):
        raise ValueError('missing archive entry')
    empty = json.loads((evidence/'empty-files.json').read_text(encoding='utf-8'))
    if empty['size'] != 0 or empty['sha256'] != sha(b'') or len(empty['files']) != manifest['zero_byte_file_count']:
        raise ValueError('changed zero-byte file index')
    reconstructed = set()
    for name in empty['files']:
        path = Path(name)
        if path.is_absolute() or '..' in path.parts or name in seen or name in reconstructed:
            raise ValueError('unsafe/duplicate zero-byte path')
        target = evidence/path; target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(b'')
        if target.stat().st_size != 0 or sha(target.read_bytes()) != empty['sha256']:
            raise ValueError('zero-byte restoration differs')
        reconstructed.add(name)
    original_count = len(seen)-1+len(reconstructed)
    if original_count != manifest['original_file_count']:
        raise ValueError('original file count differs')
    script = evidence/'sources/interim'/HERE.name/'measure_nc.py'
    nc = out/'nc-measurements.json'
    subprocess.run([sys.executable, '-I', '-X', 'utf8', str(script), '--source', str(evidence/'nc/paths'), '--output', str(nc)], check=True)
    if nc.read_bytes() != (HERE/'nc-measurements.json').read_bytes():
        raise ValueError('NC measurements differ')
    code = """import json,sys
from pathlib import Path
source,evidence=map(Path,sys.argv[1:3]);sys.path.insert(0,str(source/'tools'))
import audit_wet_joint_bulk_history as bulk
for folder in sorted((evidence/'joint').glob('g*-l8')):
    actual=bulk.audit(folder/'accepted',json.loads((folder/'result.json').read_text()),bulk.history.inputs())
    assert actual==json.loads((folder/'audit.json').read_text()),'joint audit differs'
comparison=bulk.history.compare_accepted(evidence/'joint/g16-l8/accepted',evidence/'joint/g32-l8/accepted',ratio=2)
assert comparison==json.loads((evidence/'joint/measurements.json').read_text())['comparison'],'joint sensitivity differs'
"""
    subprocess.run([sys.executable, '-I', '-X', 'utf8', '-c', code, str(evidence/'sources'), str(evidence)], check=True)
    result = dict(archive_verified=True, verified_entries=original_count, stored_archive_entries=len(seen),
                  reconstructed_zero_byte_files=len(reconstructed), verified_parts=len(manifest['parts']),
                  nc_measurements_byte_identical=True, joint_audits_and_sensitivity_exact=True,
                  sealed_sources_in_isolated_utf8_subprocess=True, numerical_acceptance=False)
    (out/'verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    print(json.dumps(run(p.parse_args().output.resolve()), indent=2))
