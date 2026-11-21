"""Package full accepted fields and solved systems, excluding material trial IO."""
import argparse
import hashlib
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def seal(evidence, destination):
    archive = evidence.parent/(evidence.name+'.tar.xz')
    entries = {p.relative_to(evidence).as_posix(): sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()}
    with tarfile.open(archive, 'w:xz', preset=3) as z:
        for name in entries:
            z.add(evidence/name, arcname=name, recursive=False)
    parts = []
    with archive.open('rb') as f:
        while block := f.read(16*1024*1024):
            name = f'evidence.tar.xz.part{len(parts)+1:02d}'
            with (destination/name).open('xb') as target:
                target.write(block)
            parts.append(dict(name=name, bytes=len(block), sha256=hashlib.sha256(block).hexdigest()))
    manifest = dict(parts=parts, archive_sha256=sha(archive), archive_bytes=archive.stat().st_size,
                    entries=entries, entry_count=len(entries), binaries_included=False,
                    coverage='all accepted full fields, contacts, actual bordered systems and Newton event records; native per-call material request/result files are hashed only, not bundled')
    (destination/'archive-sha256.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--evidence', type=Path, required=True); p.add_argument('--destination', type=Path, default=HERE)
    a = p.parse_args(); seal(a.evidence.resolve(), a.destination.resolve())
