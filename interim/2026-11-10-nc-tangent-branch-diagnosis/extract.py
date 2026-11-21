"""Verify the previous complete archive; retain only inputs to this diagnosis."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / 'interim/2026-11-10-nc-bounded-refinement'


def selected(name):
    p = Path(name)
    if name.startswith('prepared/origin/') or name.startswith('prepared/observer/'):
        return True
    if name.startswith('prepared/reference400/'):
        return True
    return (p.parts[0].startswith(('native-', 'independent-'))
            and not p.parts[0].endswith('-resumed')
            and p.name in {'controls.json', 'report.json', 'restart.json',
                           'samples.jsonl.gz', 'attempts.jsonl.gz',
                           'terminal-linear-systems.jsonl.f64.gz'})


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((PRIOR / 'archive-sha256.json').read_text())
    spec = importlib.util.spec_from_file_location('bounded_archive_reader', PRIOR / 'remeasure.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    total = hashlib.sha256()
    size = 0
    for part in manifest['parts']:
        checksum = hashlib.sha256()
        count = 0
        with (PRIOR / part['name']).open('rb') as stream:
            for raw in iter(lambda: stream.read(1024 * 1024), b''):
                checksum.update(raw)
                total.update(raw)
                count += len(raw)
                size += len(raw)
        if checksum.hexdigest() != part['sha256'] or count != part['size']:
            raise ValueError('archive part changed: ' + part['name'])
    if total.hexdigest() != manifest['archive_sha256'] or size != manifest['archive_size']:
        raise ValueError('complete archive changed')
    seen = set()
    retained = {}
    reader = module.PartsReader([PRIOR / p['name'] for p in manifest['parts']])
    try:
        with tarfile.open(fileobj=reader, mode='r|gz') as archive:
            for entry in archive:
                p = Path(entry.name)
                if (not entry.isfile() or p.is_absolute() or '..' in p.parts
                        or entry.name in seen or entry.name not in manifest['entries']):
                    raise ValueError('unsafe or unlisted archive member')
                expected = manifest['entries'][entry.name]
                if entry.size != expected['size']:
                    raise ValueError('archive member size changed')
                keep = selected(entry.name)
                target = out / entry.name
                if keep:
                    target.parent.mkdir(parents=True, exist_ok=True)
                checksum = hashlib.sha256()
                with archive.extractfile(entry) as stream:
                    destination = target.open('xb') if keep else None
                    try:
                        for raw in iter(lambda: stream.read(1024 * 1024), b''):
                            checksum.update(raw)
                            if destination is not None:
                                destination.write(raw)
                    finally:
                        if destination is not None:
                            destination.close()
                if checksum.hexdigest() != expected['sha256']:
                    raise ValueError('archive member changed: ' + entry.name)
                seen.add(entry.name)
                if keep:
                    retained[entry.name] = expected
    finally:
        reader.close()
    if seen != set(manifest['entries']):
        raise ValueError('archive inventory changed')
    result = dict(archive_sha256=total.hexdigest(), archive_size=size,
                  verified_entries=len(seen), retained_entries=retained,
                  starting_commit='113d83b8a37b3b68ce46eb00b6664f042334eae0',
                  numerical_acceptance=False)
    (out / 'extraction.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.output.resolve())
    print(json.dumps({k: v for k, v in result.items() if k != 'retained_entries'}, indent=2))
