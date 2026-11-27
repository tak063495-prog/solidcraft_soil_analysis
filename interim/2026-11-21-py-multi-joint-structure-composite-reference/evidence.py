"""Package/verify unchanged raw histories and the independent implementation."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
from reference import dump


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pack(directories, archive, manifest):
    files = {}
    for directory in map(Path, directories):
        study = json.loads((directory/'study.json').read_text(encoding='utf-8'))
        if not study['complete']:
            raise ValueError('partial study')
        for name, expected in study['raw_sha256'].items():
            path = directory/name
            key = 'raw/'+name
            if key in files or digest(path) != expected:
                raise ValueError('duplicate or modified raw history: '+name)
            files[key] = path
    here = Path(__file__).resolve().parent
    for path in sorted(here.glob('*.py')):
        files['sources/'+path.name] = path
    files['sources/requirements.txt'] = here/'requirements.txt'
    entries = {name: dict(sha256=digest(path), bytes=path.stat().st_size) for name, path in files.items()}
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for name, path in sorted(files.items()):
            z.write(path, name)
    result = dict(schema='dry-multi-joint-evidence/1', accepted=None, archive_sha256=digest(archive),
                  archive_bytes=Path(archive).stat().st_size, entries=entries,
                  numerical_acceptance=False, production_solver_or_solution_used=False)
    dump(manifest, result)
    return result


def verify(archive, manifest):
    value = json.loads(Path(manifest).read_text(encoding='utf-8'))
    if digest(archive) != value['archive_sha256']:
        raise ValueError('archive SHA differs')
    with zipfile.ZipFile(archive) as z:
        if len(z.namelist()) != len(set(z.namelist())) or set(z.namelist()) != set(value['entries']):
            raise ValueError('archive entry inventory differs')
        for name, expected in value['entries'].items():
            if Path(name).is_absolute() or '..' in Path(name).parts:
                raise ValueError('unsafe archive path')
            raw = z.read(name)
            if len(raw) != expected['bytes'] or hashlib.sha256(raw).hexdigest() != expected['sha256']:
                raise ValueError('raw entry differs: '+name)
    return dict(archive_verified=True, entries_verified=len(value['entries']),
                accepted=None, numerical_acceptance=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('pack', 'verify'))
    p.add_argument('--studies', nargs='+', type=Path)
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    if a.action == 'pack':
        if not a.studies:
            p.error('--studies required for pack')
        result = pack(a.studies, a.archive, a.manifest)
    else:
        result = verify(a.archive, a.manifest)
    if a.output:
        dump(a.output, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'entries'}))


if __name__ == '__main__':
    main()
