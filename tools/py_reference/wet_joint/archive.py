"""Authenticate split evidence and extract to a new directory (no new solve)."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import zipfile
from .compare import digest, read, write


def safe_name(name):
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or '\\' in name or ':' in name or not path.parts:
        raise ValueError('unsafe archive path: '+name)
    return path


def extract(index_path, output):
    index_path = Path(index_path); output = Path(output)
    index = read(index_path); output.mkdir(parents=True, exist_ok=False)
    combined = output/'evidence.zip'
    try:
        with combined.open('xb') as dest:
            for part in index['parts']:
                if Path(part['file']).name != part['file']:
                    raise ValueError('unsafe evidence part')
                path = index_path.parent/part['file']
                if path.stat().st_size != part['bytes'] or digest(path) != part['sha256']:
                    raise ValueError('evidence part SHA mismatch')
                with path.open('rb') as src:
                    shutil.copyfileobj(src, dest)
        if digest(combined) != index['archive_sha256']:
            raise ValueError('whole archive SHA mismatch')
        expected = {v['file']: v for v in index['entries']}
        if len(expected) != len(index['entries']):
            raise ValueError('duplicate index entry')
        with zipfile.ZipFile(combined) as archive:
            if len(archive.namelist()) != len(expected) or set(archive.namelist()) != set(expected):
                raise ValueError('archive inventory mismatch')
            # Verify every byte before extracting any member.
            for name in archive.namelist():
                safe_name(name); record = expected[name]; h = hashlib.sha256(); count = 0
                with archive.open(name) as src:
                    while block := src.read(1024*1024):
                        h.update(block); count += len(block)
                if h.hexdigest() != record['sha256'] or count != record['bytes']:
                    raise ValueError('archive entry SHA mismatch: '+name)
            for name in archive.namelist():
                path = output/'data'/str(safe_name(name)); path.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(name) as src, path.open('xb') as dst:
                    shutil.copyfileobj(src, dst)
        report = dict(software_integrity_passed=True, entries=len(expected), numerical_acceptance=False)
        write(output/'verification.json', report)
        combined.unlink()
        return report
    except BaseException as error:
        write(output/'failure.json', dict(completed=False, passed=False, error=repr(error)))
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--index', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    print(json.dumps(extract(a.index, a.output)))


if __name__ == '__main__':
    main()
