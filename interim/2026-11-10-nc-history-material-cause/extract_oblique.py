"""Verify the complete existing oblique archive; retain only stress-path inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

HERE = Path(__file__).resolve().parent
PRIOR = HERE.with_name('2026-11-10-nc-oblique-and-joint-history')


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((PRIOR/'archive-sha256.json').read_text(encoding='utf-8'))
    archive = out/'joined.zip'; digest = hashlib.sha256(); size = 0
    with archive.open('xb') as stream:
        for item in manifest['parts']:
            raw = (PRIOR/item['name']).read_bytes()
            if len(raw) != item['size'] or sha(raw) != item['sha256']:
                raise ValueError('prior archive part mismatch')
            stream.write(raw); digest.update(raw); size += len(raw)
    if digest.hexdigest() != manifest['archive_sha256'] or size != manifest['archive_size']:
        raise ValueError('prior archive mismatch')
    kept = {}; seen = set()
    with zipfile.ZipFile(archive) as z:
        for item in z.infolist():
            name = item.filename; path = Path(name)
            if path.is_absolute() or '..' in path.parts or item.is_dir() or name in seen:
                raise ValueError('unsafe prior member')
            expected = manifest['entries'].get(name); raw = z.read(item)
            if expected is None or len(raw) != expected['size'] or sha(raw) != expected['sha256']:
                raise ValueError('prior member mismatch')
            seen.add(name)
            keep = (name.startswith('nc/paths/') and 'a80-b20' in name and '-resumed/' not in name
                    and path.name in {'model.json', 'samples.jsonl.gz', 'report.json', 'attempts.jsonl.gz'})
            if keep:
                target = out/path; target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(raw)
                kept[name] = expected
    if seen != set(manifest['entries']):
        raise ValueError('prior member inventory mismatch')
    archive.unlink()
    report = dict(archive_sha256=digest.hexdigest(), verified_parts=len(manifest['parts']),
                  verified_stored_entries=len(seen), retained=kept, numerical_acceptance=False,
                  scope='All stored members verified; unrelated zero-byte process logs not reconstructed.')
    (out/'extraction.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    r = run(p.parse_args().output.resolve())
    print(json.dumps({k: v for k, v in r.items() if k != 'retained'}, indent=2))
