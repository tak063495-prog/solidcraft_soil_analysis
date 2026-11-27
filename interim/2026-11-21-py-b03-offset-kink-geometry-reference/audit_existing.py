"""Read-only audit: existing index OIDs/modes and available working bytes vs base.

Sparse omitted paths are reported, never claimed as re-read or revalidated.
"""
import argparse
import hashlib
import json
from pathlib import Path
import stat
import subprocess

BASE_SHA = 'f377e8eda41777c84c3bb849af556bf61683b4cb'


def audit(root):
    root = Path(root).resolve()
    raw = subprocess.check_output(['git', 'ls-tree', '-r', '-z', BASE_SHA], cwd=root)
    index_raw = subprocess.check_output(['git', 'ls-files', '--stage', '-z'], cwd=root)
    flags_raw = subprocess.check_output(['git', 'ls-files', '-v', '-z'], cwd=root)
    sparse = {r[2:].decode() for r in flags_raw.split(b'\0') if r.startswith(b'S ')}
    index = {}
    for record in index_raw.split(b'\0'):
        if record:
            header, name = record.split(b'\t', 1)
            mode, oid, stage = header.decode().split()
            index[name.decode()] = (mode, oid, stage)
    mismatches, omitted, count, checked, total = [], [], 0, 0, 0
    for record in raw.split(b'\0'):
        if not record:
            continue
        header, name = record.split(b'\t', 1)
        mode, kind, oid = header.decode().split()
        name = name.decode()
        path = root/name
        count += 1
        if index.get(name) != (mode, oid, '0'):
            mismatches.append({'path': name, 'reason': 'existing index OID/mode differs from base'})
        if not path.is_file():
            if name in sparse:
                omitted.append({'path': name, 'mode': mode, 'git_oid': oid})
            else:
                mismatches.append({'path': name, 'reason': 'missing non-sparse tracked file'})
            continue
        if kind != 'blob':
            mismatches.append({'path': name, 'reason': 'not an ordinary tracked blob'})
            continue
        size = path.stat().st_size
        checked += 1
        total += size
        digest = hashlib.sha1(f'blob {size}\0'.encode())
        with path.open('rb') as f:
            for chunk in iter(lambda: f.read(1024*1024), b''):
                digest.update(chunk)
        is_exec = bool(path.stat().st_mode & stat.S_IXUSR)
        if digest.hexdigest() != oid or is_exec != (mode == '100755'):
            mismatches.append({'path': name, 'expected_oid': oid, 'actual_oid': digest.hexdigest(),
                               'expected_mode': mode, 'actual_executable': is_exec})
    return {'base_sha': BASE_SHA, 'existing_index_paths_checked': count,
            'working_copy_files_checked': checked, 'working_copy_bytes_checked': total,
            'tree_listing_sha256': hashlib.sha256(raw).hexdigest(), 'mismatches': mismatches,
            'sparse_omitted_count': len(omitted), 'sparse_omitted_paths': omitted,
            'existing_tree_preserved_and_available_bytes_identical': not mismatches,
            'omitted_file_contents_revalidated': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository_root')
    parser.add_argument('output')
    args = parser.parse_args()
    result = audit(args.repository_root)
    Path(args.output).write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k != 'sparse_omitted_paths'}, ensure_ascii=False))
    raise SystemExit(0 if result['existing_tree_preserved_and_available_bytes_identical'] else 1)
