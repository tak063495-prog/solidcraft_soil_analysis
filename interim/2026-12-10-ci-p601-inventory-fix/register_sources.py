#!/usr/bin/env python3
"""Print a pinned TOML source registration from committed, unchanged JSON bytes.

Append the output to inventory.toml (keep its single schema line), then update
SUPPLEMENTAL_SHA256 in audit_status.py. Never regenerate old source records.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess


def status_values(value, pointer=''):
    if isinstance(value, dict):
        for key, child in value.items():
            token = key.replace('~', '~0').replace('/', '~1')
            yield from status_values(child, pointer+'/'+token)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from status_values(child, pointer+'/'+str(index))
    elif isinstance(value, bool) or value is None:
        yield pointer, value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--commit', required=True)
    parser.add_argument('--id', required=True)
    parser.add_argument('path')
    args = parser.parse_args()
    if not re.fullmatch(r'[a-f0-9]{40}', args.commit):
        parser.error('--commit requires a full immutable SHA')
    root = Path(__file__).resolve().parents[2]
    path = Path(args.path)
    if path.is_absolute() or '..' in path.parts or not args.path.startswith('interim/'):
        parser.error('path must be a repository-relative interim JSON')
    frozen = subprocess.check_output(['git', '-C', str(root), 'show', args.commit+':'+args.path])
    if (root/path).read_bytes() != frozen:
        parser.error('working bytes differ from the specified commit')
    proofs = list(status_values(json.loads(frozen)))
    if not proofs:
        parser.error('JSON must contain typed boolean/null status evidence')
    source = dict(id=args.id, path=args.path, commit=args.commit, kind='json',
                  sha256=hashlib.sha256(frozen).hexdigest(),
                  git_blob_sha1=hashlib.sha1(b'blob '+str(len(frozen)).encode()+b'\0'+frozen).hexdigest())
    print('[[sources]]')
    for key, value in source.items():
        print(key+' = '+json.dumps(value, ensure_ascii=False))
    for pointer, value in proofs:
        print('[[sources.status_keys]]')
        print('pointer = '+json.dumps(pointer, ensure_ascii=False))
        print('value_json = '+json.dumps(json.dumps(value, ensure_ascii=False), ensure_ascii=False))
    print()


if __name__ == '__main__':
    main()
