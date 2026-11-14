"""Verify the prior archive in full, retaining only the independent original seed.

No material or equilibrium solve; archive files and old records are immutable.
"""
import argparse
import hashlib
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / '2026-11-09-nc-pseudo-arclength'


def prepare(out):
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((PRIOR/'archive-sha256.json').read_text())
    archive = out/'joined.tar.xz'
    h = hashlib.sha256()
    with archive.open('xb') as f:
        for part in manifest['parts']:
            data = (PRIOR/part['name']).read_bytes()
            if len(data) != part['bytes'] or hashlib.sha256(data).hexdigest() != part['sha256']:
                raise ValueError('archive part differs')
            f.write(data)
            h.update(data)
    if h.hexdigest() != manifest['archive_sha256']:
        raise ValueError('archive differs')
    seen = set()
    retained = {}
    with tarfile.open(archive, 'r|xz') as z:
        for m in z:
            name = Path(m.name)
            if not m.isfile() or name.is_absolute() or '..' in name.parts or m.name in seen:
                raise ValueError('unsafe/duplicate entry')
            data = z.extractfile(m).read()
            actual = hashlib.sha256(data).hexdigest()
            if manifest['entries'].get(m.name) != actual:
                raise ValueError('entry differs: '+m.name)
            seen.add(m.name)
            if m.name == 'independent-arc-history400/samples.json':
                dest = out/m.name
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
                retained[m.name] = actual
    if seen != set(manifest['entries']):
        raise ValueError('archive inventory differs')
    report = dict(passed=False, solver_executed=False, material_map_executed=False,
                  archive_sha256=h.hexdigest(), verified_entries=len(seen), retained=retained,
                  note='Full byte verification; not a rerun of prior numerical measurements')
    (out/'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    print(json.dumps(prepare(p.parse_args().output.resolve()), indent=2))
