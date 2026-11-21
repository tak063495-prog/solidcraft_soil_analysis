"""Bounded SHA verification and split evidence packaging; no numerical reruns."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def verify(folder, extract=None):
    folder = Path(folder); m = json.loads((folder/'archive-index.json').read_text(encoding='utf-8'))
    h = hashlib.sha256(); size = 0
    if extract:
        extract.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryFile() as stream:
        for part in m['parts']:
            p = folder/part['name']
            if p.stat().st_size != part['size'] or sha(p) != part['sha256']:
                raise ValueError('evidence part SHA/size mismatch: '+str(p))
            with p.open('rb') as f:
                for b in iter(lambda: f.read(1024*1024), b''):
                    stream.write(b); h.update(b); size += len(b)
        if size != m['archive_bytes'] or h.hexdigest() != m['archive_sha256']:
            raise ValueError('combined archive differs')
        stream.seek(0)
        with zipfile.ZipFile(stream) as z:
            names = z.namelist()
            if len(names) != len(set(names)) or set(names) != set(m['entries']):
                raise ValueError('archive member inventory differs')
            for name in names:
                p = Path(name)
                if p.is_absolute() or '..' in p.parts:
                    raise ValueError('unsafe archive member')
                raw = z.read(name); e = m['entries'][name]
                if len(raw) != e['size'] or hashlib.sha256(raw).hexdigest() != e['sha256']:
                    raise ValueError('member SHA/size differs: '+name)
                if extract:
                    target = extract/p; target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(raw)
    return dict(parts=len(m['parts']), entries=len(m['entries']), all_sha256_verified=True,
                numerical_reruns=0)


def seal(destination, folders):
    destination = Path(destination); entries = {}
    with tempfile.TemporaryFile() as stream:
        with zipfile.ZipFile(stream, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            for folder in folders:
                for p in sorted(folder.rglob('*')):
                    if not p.is_file() or '__pycache__' in p.parts:
                        continue
                    name = folder.name+'/'+p.relative_to(folder).as_posix()
                    if name in entries:
                        raise ValueError('duplicate member')
                    z.write(p, name); entries[name] = dict(size=p.stat().st_size, sha256=sha(p))
        stream.seek(0); parts = []; h = hashlib.sha256(); total = 0
        while raw := stream.read(16*1024*1024):
            name = 'evidence.zip.part'+str(len(parts)).zfill(3); p = destination/name
            with p.open('xb') as f:
                f.write(raw)
            digest = hashlib.sha256(raw).hexdigest(); h.update(raw); total += len(raw)
            parts.append(dict(name=name, size=len(raw), sha256=digest))
    (destination/'archive-index.json').write_text(json.dumps(dict(
        schema='c04follow-evidence/1', archive_sha256=h.hexdigest(), archive_bytes=total,
        parts=parts, entries=entries), indent=2)+'\n', encoding='utf-8')
    return verify(destination)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verify', type=Path)
    p.add_argument('--extract', type=Path)
    p.add_argument('--seal', type=Path)
    p.add_argument('--folder', action='append', type=Path, default=[])
    a = p.parse_args()
    print(verify(a.verify, a.extract) if a.verify else seal(a.seal, a.folder))
