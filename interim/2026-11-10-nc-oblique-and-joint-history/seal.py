"""Package new evidence and the exact read-only audit source snapshot."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def run(args):
    archive = args.output/'new-evidence.zip'
    args.output.mkdir(parents=True, exist_ok=False)
    entries = {}
    empty = set()
    def add(z, name, raw):
        if name in entries or name in empty:
            raise ValueError('duplicate evidence member')
        if not raw:
            empty.add(name)
            return
        entries[name] = dict(size=len(raw), sha256=sha(raw))
        z.writestr(name, raw)
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for source, prefix in ((args.nc, 'nc'), (args.joint, 'joint'), (args.prefix, 'fresh-prefix')):
            for path in sorted(source.rglob('*')):
                if path.is_file():
                    add(z, prefix+'/'+path.relative_to(source).as_posix(), path.read_bytes())
        if args.failed is not None:
            for path in sorted(args.failed.rglob('*')):
                if path.is_file():
                    add(z, 'recording-failure/'+path.relative_to(args.failed).as_posix(), path.read_bytes())
        for path in args.provenance:
            add(z, 'provenance/'+path.name, path.read_bytes())
        snapshot = json.loads(args.sources.read_text(encoding='utf-8'))
        for name, expected in snapshot.items():
            raw = (ROOT/name).read_bytes()
            if sha(raw) != expected:
                raise ValueError('tested source changed: '+name)
            add(z, 'sources/'+name, raw)
        for path in sorted((ROOT/'interim').rglob('*')):
            if path.is_file() and path.suffix in ('.py', '.rs') and 'work' not in path.parts and '__pycache__' not in path.parts:
                add(z, 'sources/'+path.relative_to(ROOT).as_posix(), path.read_bytes())
        add(z, 'empty-files.json', (json.dumps(dict(size=0, sha256=sha(b''), files=sorted(empty)), indent=2)+'\n').encode())
    raw = archive.read_bytes(); parts = []
    for offset in range(0, len(raw), 2*1024*1024):
        name = f'evidence.zip.part{len(parts):03d}'
        part = raw[offset:offset+2*1024*1024]
        (HERE/name).write_bytes(part)
        parts.append(dict(name=name, size=len(part), sha256=sha(part)))
    manifest = dict(archive_sha256=sha(raw), archive_size=len(raw), entries=entries, parts=parts,
                    original_file_count=len(entries)-1+len(empty), zero_byte_file_count=len(empty),
                    numerical_acceptance=False)
    (HERE/'archive-sha256.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('nc', 'joint', 'prefix', 'sources', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--provenance', type=Path, nargs='+', required=True)
    p.add_argument('--failed', type=Path)
    run(p.parse_args())
