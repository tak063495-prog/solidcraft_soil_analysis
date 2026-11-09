"""Independently extract/verify/re-measure, without a temporary joined archive.

This storage helper is versioned beside the evidence archive; the seven archived
numerical/audit helper snapshots are also verified. No source/cache is deleted.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile

from measure import measure

HERE = Path(__file__).resolve().parent


class Parts:
    def __init__(self, paths):
        self.paths = iter(paths)
        self.current = None

    def read(self, n):
        if n < 0:
            raise ValueError('bounded stream read required')
        result = bytearray()
        while len(result) < n:
            if self.current is None:
                try:
                    self.current = next(self.paths).open('rb')
                except StopIteration:
                    break
            block = self.current.read(n-len(result))
            if block:
                result.extend(block)
            else:
                self.current.close()
                self.current = None
        return bytes(result)


def remeasure(out):
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((HERE/'archive-sha256.json').read_text())
    archive_hash = hashlib.sha256()
    for part in manifest['parts']:
        h = hashlib.sha256()
        count = 0
        with (HERE/part['name']).open('rb') as f:
            while block := f.read(1024*1024):
                archive_hash.update(block); h.update(block); count += len(block)
        if count != part['bytes'] or h.hexdigest() != part['sha256']:
            raise ValueError('archive part differs: '+part['name'])
    if archive_hash.hexdigest() != manifest['archive_sha256']:
        raise ValueError('whole archive differs')
    evidence = out/'evidence'
    evidence.mkdir()
    seen = set()
    with tarfile.open(fileobj=Parts([HERE/p['name'] for p in manifest['parts']]), mode='r|xz') as z:
        for member in z:
            name = PurePosixPath(member.name)
            if not member.isfile() or name.is_absolute() or '..' in name.parts or ':' in name.parts[0] or member.name in seen:
                raise ValueError('unsafe or duplicate archive entry')
            dest = evidence/member.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            h = hashlib.sha256()
            with z.extractfile(member) as src, dest.open('xb') as target:
                while block := src.read(1024*1024):
                    h.update(block); target.write(block)
            if h.hexdigest() != manifest['entries'].get(member.name):
                raise ValueError('entry differs: '+member.name)
            seen.add(member.name)
    if seen != set(manifest['entries']):
        raise ValueError('archive inventory differs')
    for name, expected in json.loads((evidence/'diagnostic-source-sha256.json').read_text()).items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != expected:
            raise ValueError('archived source differs: '+name)
    result = measure(evidence, out/'recomputed')
    saved = HERE/'measurements.json'
    computed = out/'recomputed/measurements.json'
    record = dict(passed=False, solver_executed=False, material_map_executed=False,
        archive_entries_verified=len(seen), archive_sha256=archive_hash.hexdigest(),
        exact_saved_measurements=result == json.loads(saved.read_text()),
        exact_saved_measurement_bytes=computed.read_bytes() == saved.read_bytes(),
        independent_extraction=True, hardlinks_used=False, temporary_joined_archive_created=False,
        streaming_driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        streaming_driver_location='interim/2026-11-09-nc-continuation-ci/stream_remeasure.py',
        streaming_driver_bundled_in_archive=False)
    (out/'remeasurement.json').write_text(json.dumps(record, indent=2)+'\n')
    if not record['exact_saved_measurements'] or not record['exact_saved_measurement_bytes']:
        raise ValueError('remeasurement differs; retain both')
    return record


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(remeasure(a.output.resolve()), indent=2))
