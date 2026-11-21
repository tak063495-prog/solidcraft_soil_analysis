"""Stream baseline byte checks without lazy-fetching historical archive blobs."""
import argparse
import hashlib
import os
from pathlib import Path
import subprocess
from a03follow_binding import ROOT
from a03follow_run import write

BASE = '89e9dfc26cbe444199179d7faa1edb1437849309'
ALLOWED = {'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml',
           'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'}


def audit(output, allow_registration=False):
    env = dict(os.environ, GIT_NO_LAZY_FETCH='1')
    files, absent, changed = {}, [], []
    for line in subprocess.check_output(['git', 'ls-tree', '-r', BASE], cwd=ROOT, env=env, text=True).splitlines():
        meta, name = line.split('\t')
        mode, kind, oid = meta.split()
        if kind != 'blob':
            continue
        path = ROOT/name
        if not path.exists():
            absent.append(name)
            continue
        size = path.stat().st_size
        blob, digest = hashlib.sha1(('blob '+str(size)+'\0').encode()), hashlib.sha256()
        with path.open('rb') as f:
            for raw in iter(lambda: f.read(1024*1024), b''):
                blob.update(raw)
                digest.update(raw)
        actual = blob.hexdigest()
        files[name] = dict(size=size, git_blob_sha1=actual, sha256=digest.hexdigest(), baseline_blob=oid)
        if actual != oid:
            changed.append(name)
    assert all(p.startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part') for p in absent)
    assert set(changed) <= (ALLOWED if allow_registration else set()), changed
    write(output, dict(base=BASE, files=files, available=len(files), baseline_byte_changes=changed,
                       missing_only_existing_ci_postfold_exclusion=absent, acceptance_values_or_archives_changed=False,
                       permitted_changes_are_only_inventory_append_and_auditor_pin=allow_registration))
    print('baseline available', len(files), 'changes', changed, 'excluded', len(absent), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--allow-registration', action='store_true')
    a = p.parse_args()
    audit(a.output, a.allow_registration)
