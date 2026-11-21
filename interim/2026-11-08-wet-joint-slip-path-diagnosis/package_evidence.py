"""Archive the recorded small-case evidence in GitHub-sized immutable shards."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
from slip_study import HERE, write, sha


def package(raw):
    entries = {}
    archives = []
    archive = None
    size = 0
    try:
        for path in sorted(raw.rglob('*')):
            if not path.is_file() or '__pycache__' in path.parts:
                continue
            data = path.read_bytes()
            if archive is None or size+len(data) > 23_000_000:
                if archive:
                    archive.close()
                target = HERE/f'evidence-{len(archives)+1:03d}.zip'
                assert not target.exists()
                archive = zipfile.ZipFile(target, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6)
                archives.append(target)
                size = 0
            name = 'raw/'+path.relative_to(raw).as_posix()
            archive.writestr(name, data)
            entries[name] = hashlib.sha256(data).hexdigest()
            size += len(data)
    finally:
        if archive:
            archive.close()
    write(HERE/'archive-sha256.json', dict(schema='wet-joint-slip-evidence/1', entries=entries,
          archives=[dict(file=p.name, sha256=sha(p), bytes=p.stat().st_size) for p in archives]))
    return dict(archives=len(archives), entries=len(entries), bytes=sum(p.stat().st_size for p in archives))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--raw',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(package(args.raw.resolve()),indent=2))
