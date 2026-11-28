"""Verify and extract every evidence shard to a new directory; no solving."""
import argparse
import json
from pathlib import Path, PurePosixPath
import zipfile
import study


def restore(output):
    manifest = json.loads((study.HERE/'archive-sha256.json').read_text())
    # Check all archive hashes before creating any extracted evidence.
    for archive in manifest['archives']:
        assert study.sha((study.HERE/archive['file']).read_bytes()) == archive['sha256'], archive['file']
    output.mkdir(parents=True, exist_ok=False)
    found = set()
    for archive in manifest['archives']:
        with zipfile.ZipFile(study.HERE/archive['file']) as stream:
            for name in stream.namelist():
                path = PurePosixPath(name)
                assert not path.is_absolute() and '..' not in path.parts and '\\' not in name
                assert name in manifest['entries'] and name not in found, ('unexpected entry', name)
                raw = stream.read(name)
                assert study.sha(raw) == manifest['entries'][name], ('changed entry', name)
                target = output.joinpath(*path.parts); target.parent.mkdir(parents=True, exist_ok=True)
                with target.open('xb') as file:
                    file.write(raw)
                found.add(name)
    assert found == set(manifest['entries']), 'missing evidence'
    return dict(all_archive_hashes_verified=True, verified_entries=len(found), solver_rerun=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(restore(args.output.resolve())))
