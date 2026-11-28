"""Authenticate the new diagnostic archive; no global solves or old rewrites."""
import hashlib
import json
from pathlib import Path
import zipfile
import shutil
import tempfile

HERE=Path(__file__).resolve().parent


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def verify():
    index=json.loads((HERE/'archive-index.json').read_text(encoding='utf-8'))
    acceptance=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
    if acceptance['passed'] is not False or acceptance['archive_index_sha256']!=digest(HERE/'archive-index.json'):
        raise ValueError('C-01 status/index changed')
    with tempfile.TemporaryDirectory() as tmp:
        archive=Path(tmp)/'evidence.zip'
        with archive.open('xb') as out:
            for part in index['parts']:
                path=HERE/part['file']
                if path.stat().st_size!=part['bytes'] or digest(path)!=part['sha256']:raise ValueError('diagnostic part integrity')
                with path.open('rb') as f:shutil.copyfileobj(f,out)
        if digest(archive)!=index['sha256'] or archive.stat().st_size!=index['bytes']:
            raise ValueError('diagnostic archive integrity')
        expected={r['path']:r for r in index['entries']}
        with zipfile.ZipFile(archive) as z:
            if set(z.namelist())!=set(expected) or len(z.namelist())!=len(expected):raise ValueError('diagnostic ZIP inventory')
            for name in z.namelist():
                p=Path(name)
                if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe ZIP entry')
                h=hashlib.sha256();count=0
                with z.open(name) as f:
                    for b in iter(lambda:f.read(1<<20),b''):h.update(b);count+=len(b)
                if h.hexdigest()!=expected[name]['sha256'] or count!=expected[name]['bytes']:raise ValueError('diagnostic entry integrity: '+name)
            for path in sorted(HERE.glob('c01_*.py')):
                name='source/'+path.name
                if name not in expected or digest(path)!=expected[name]['sha256']:raise ValueError('sealed diagnostic source differs: '+path.name)
    return dict(archive_verified=True,new_global_equilibrium_solved=False,c01_passed=False)

if __name__=='__main__':print(json.dumps(verify()))
