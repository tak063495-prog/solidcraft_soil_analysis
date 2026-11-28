"""Stream-hash every present baseline file without downloading omitted blobs."""
import argparse
import hashlib
from pathlib import Path
import subprocess
from b01_cases import ROOT, HERE, read, write

BASE = '3a20ac2974c46f60dde56fae878826962a41d3be'
DOCS = {'docs/P6_STATUS.ja.md', 'docs/STRUCTURE_WATER_JOINT_REMAINDERS.ja.md',
        'interim/2026-11-16-issue-triage/todo.ja.md'}


def run(study, output):
    rows, absent, prefixes = [], [], []
    tree = subprocess.check_output(['git', 'ls-tree', '-r', BASE], cwd=ROOT, text=True)
    for line in tree.splitlines():
        meta, name = line.split('\t')
        mode, kind, blob = meta.split()
        path = ROOT/name
        if not path.is_file():
            absent.append(dict(path=name, mode=mode, blob=blob))
            continue
        actual_mode = '100755' if path.stat().st_mode & 0o111 else '100644'
        if name in DOCS:
            original = subprocess.check_output(['git', 'show', BASE+':'+name], cwd=ROOT)
            current = path.read_bytes()
            prefixes.append(dict(path=name, original_bytes=len(original), original_sha256=hashlib.sha256(original).hexdigest(),
                current_bytes=len(current), current_sha256=hashlib.sha256(current).hexdigest(),
                all_original_bytes_preserved_as_prefix=current.startswith(original), mode_unchanged=mode == actual_mode))
            assert current.startswith(original) and mode == actual_mode
        else:
            size = path.stat().st_size
            digest = hashlib.sha1(f'blob {size}\0'.encode())
            with path.open('rb') as f:
                for block in iter(lambda: f.read(1024*1024), b''):
                    digest.update(block)
            assert digest.hexdigest() == blob and actual_mode == mode, name
            rows.append(dict(path=name, bytes=size, blob=blob, mode=mode))
    before, after = [read(study/f'source-binding-{name}.json') for name in ('before', 'after')]
    assert before == after
    for name, digest in before.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    write(output, dict(passed=True, baseline=BASE, baseline_paths=len(tree.splitlines()),
        present_unchanged_files=len(rows), present_unchanged_bytes=sum(r['bytes'] for r in rows),
        verified_present_files=rows, prefix_only_documents=prefixes,
        omitted_count=len(absent), omitted_original_mode_blob=absent,
        omitted_contents_not_reaudited=True, bound_fresh_sources=len(before),
        fresh_sources_equal_before_after_and_at_preservation=True,
        old_hash_snapshot_failure_and_NC_records_unchanged=True))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('study', type=Path)
    p.add_argument('output', type=Path, default=HERE/'preservation.json')
    a = p.parse_args()
    run(a.study, a.output)
