"""Verify archive bytes and remeasure all saved 3D evidence without solving."""
import argparse
import hashlib
import json
import tarfile
from pathlib import Path

from measure_further import measure

HERE = Path(__file__).resolve().parent


def remeasure(out):
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((HERE/'archive-sha256.json').read_text())
    archive = out/'joined.tar.xz'
    h = hashlib.sha256()
    with archive.open('xb') as f:
        for item in manifest['parts']:
            b = (HERE/item['name']).read_bytes()
            if len(b) != item['bytes'] or hashlib.sha256(b).hexdigest() != item['sha256']:
                raise ValueError('archive part differs')
            f.write(b); h.update(b)
    if h.hexdigest() != manifest['archive_sha256']:
        raise ValueError('joined archive differs')
    evidence = out/'evidence'; evidence.mkdir()
    with tarfile.open(archive, 'r:xz') as z:
        members = z.getmembers()
        if len(members) != len(manifest['entries']) or {m.name for m in members} != set(manifest['entries']):
            raise ValueError('archive inventory differs')
        for m in members:
            name = Path(m.name)
            if not m.isfile() or name.is_absolute() or '..' in name.parts:
                raise ValueError('unsafe archive entry')
            b = z.extractfile(m).read()
            if hashlib.sha256(b).hexdigest() != manifest['entries'][m.name]:
                raise ValueError('entry differs: '+m.name)
            p = evidence/name; p.parent.mkdir(parents=True, exist_ok=True); p.write_bytes(b)
    for name, expected in json.loads((evidence/'diagnostic-source-sha256.json').read_text()).items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != expected:
            raise ValueError('diagnostic source differs: '+name)
    result = measure(evidence, out/'recomputed')
    exact = result == json.loads((HERE/'measurements.json').read_text())
    record = dict(passed=False, solver_executed=False, material_map_executed=False,
                  archive_entries_verified=len(manifest['entries']), exact_saved_measurements=exact)
    (out/'remeasurement.json').write_text(json.dumps(record, indent=2)+'\n')
    if not exact:
        raise ValueError('remeasurement differs; retain both')
    return record


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); print(json.dumps(remeasure(a.output.resolve()), indent=2))
