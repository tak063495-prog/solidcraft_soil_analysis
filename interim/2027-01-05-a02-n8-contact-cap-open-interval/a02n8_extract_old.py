"""Copy identified old controls with part/archive/member SHA verification; no solve."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile
from a02n8_common import write, sha

BASE = '612f0f64c704f0ef0dee950a81a10f7e312a745c'


def selected(name):
    if not name.endswith(('-fields.jsonl', '/case.json', '/binding.json', '/mesh.json', '/ground.json')):
        return False
    baseline = any(name.startswith(f'formal/cross/n{n}-m4/stress_free-cycle/') for n in (4, 8))
    normal = name.startswith('formal/exchange/stress_free-cycle/')
    cap = name.startswith('retained/a02follow-contact-exchange/')
    return baseline or normal or cap


def extract(source, output):
    output.mkdir(parents=True, exist_ok=False)
    index = json.loads((source/'archive-index.json').read_text(encoding='utf-8'))
    copied = []
    digest = hashlib.sha256()
    with tempfile.TemporaryFile() as stream:
        for part in index['parts']:
            path = source/part['name']
            assert path.stat().st_size == part['size'] and sha(path) == part['sha256']
            with path.open('rb') as f:
                while raw := f.read(1024*1024):
                    stream.write(raw)
                    digest.update(raw)
        assert stream.tell() == index['archive_bytes'] and digest.hexdigest() == index['archive_sha256']
        stream.seek(0)
        with zipfile.ZipFile(stream) as archive:
            assert len(archive.namelist()) == len(set(archive.namelist()))
            for name in sorted(n for n in index['entries'] if selected(n)):
                path = Path(name)
                assert not path.is_absolute() and '..' not in path.parts
                raw = archive.read(name)
                entry = index['entries'][name]
                assert len(raw) == entry['size'] and hashlib.sha256(raw).hexdigest() == entry['sha256']
                target = output/name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(raw)
                copied.append(dict(member=name, **entry))
    (output/'source').mkdir()
    path = source/'contact-policy-exchange.json'
    (output/'source'/path.name).write_bytes(path.read_bytes())
    write(output/'input-binding.json', dict(base=BASE, source=str(source.relative_to(source.parents[1])),
        old_archive_index_sha256=sha(source/'archive-index.json'),
        old_archive_sha256=index['archive_sha256'], parts=index['parts'], members=copied,
        old_contact_exchange_summary_sha256=sha(path),
        all_selected_members_verified=True, old_files_changed=False, fresh_solves=0))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    extract(a.source.resolve(), a.output.resolve())
