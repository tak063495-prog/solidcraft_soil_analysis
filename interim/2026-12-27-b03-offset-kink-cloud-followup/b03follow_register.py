"""Append already committed JSON to the existing fixed-commit P6 inventory."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def register(commit, names):
    if not re.fullmatch('[0-9a-f]{40}', commit) or not names or len(set(names)) != len(names):
        raise ValueError('one fixed full commit and unique JSON names required')
    inventory = ROOT/'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml'
    auditor = ROOT/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'
    before = inventory.read_bytes()
    tail = ''
    for name in names:
        if name not in ('acceptance', 'measurements', 'validation'):
            raise ValueError('unexpected status document')
        path = (HERE/(name+'.json')).relative_to(ROOT).as_posix()
        raw = (ROOT/path).read_bytes()
        fixed = subprocess.check_output(['git', 'show', commit+':'+path], cwd=ROOT)
        if fixed != raw:
            raise ValueError('working bytes differ from fixed commit: '+path)
        obj = json.loads(raw)
        sid = 'b03follow_20261227_'+name
        if ('id = '+json.dumps(sid)).encode() in before or ('path = '+json.dumps(path)).encode() in before:
            raise ValueError('already registered')
        blob = hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        values = dict(id=sid, path=path, commit=commit, kind='json',
            sha256=hashlib.sha256(raw).hexdigest(), git_blob_sha1=blob)
        tail += '\n[[sources]]\n'+''.join(key+' = '+json.dumps(value)+'\n' for key, value in values.items())
        keys = ['passed', 'b03_complete', 'diagnostic_only']
        if name != 'validation':
            keys += ['general_joint_acceptance', 'production_scope_expanded']
        for key in keys:
            value = obj[key]
            if value is not None and type(value) is not bool:
                raise ValueError('status must be typed bool/null')
            if key == 'b03_complete' and value is not False:
                raise ValueError('B-03 must remain incomplete')
            tail += '[[sources.status_keys]]\npointer = '+json.dumps('/'+key)+'\nvalue_json = '+json.dumps(json.dumps(value))+'\n'
    after = before+tail.encode()
    old_auditor = auditor.read_bytes()
    new_auditor, count = re.subn(rb"SUPPLEMENTAL_SHA256 = '[0-9a-f]{64}'",
        b"SUPPLEMENTAL_SHA256 = '"+hashlib.sha256(after).hexdigest().encode()+b"'", old_auditor)
    if count != 1:
        raise ValueError('expected exactly one supplemental SHA constant')
    inventory.write_bytes(after)
    auditor.write_bytes(new_auditor)
    return dict(registered=names, commit=commit, manifest_sha256=hashlib.sha256(after).hexdigest(),
                old_inventory_is_exact_prefix=after.startswith(before), auditor_only_pinned_sha_changed=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--commit', required=True)
    p.add_argument('names', nargs='+')
    a = p.parse_args()
    print(json.dumps(register(a.commit, a.names), indent=2))
