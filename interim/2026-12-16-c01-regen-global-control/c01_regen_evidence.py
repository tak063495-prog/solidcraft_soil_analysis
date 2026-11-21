"""Lightweight authentication of sealed heavy evidence; never reruns it."""
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import zipfile

HERE=Path(__file__).resolve().parent


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def verify():
    index=json.loads((HERE/'archive-index.json').read_text(encoding='utf-8'))
    status=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
    if status['passed'] is not False or status['archive_index_sha256']!=digest(HERE/'archive-index.json'):
        raise ValueError('status or archive index changed')
    expected={e['path']:e for e in index['entries']}
    with tempfile.TemporaryDirectory() as tmp:
        archive=Path(tmp)/'evidence.zip'
        with archive.open('xb') as out:
            for p in index['parts']:
                source=HERE/p['file']
                if source.stat().st_size!=p['bytes'] or digest(source)!=p['sha256']:raise ValueError('sealed part SHA differs')
                with source.open('rb') as f:shutil.copyfileobj(f,out)
        if archive.stat().st_size!=index['bytes'] or digest(archive)!=index['sha256']:raise ValueError('sealed ZIP SHA differs')
        with zipfile.ZipFile(archive) as z:
            if len(z.namelist())!=len(expected) or set(z.namelist())!=set(expected):raise ValueError('sealed ZIP inventory differs')
            for name in sorted(expected):
                path=Path(name)
                if path.is_absolute() or '..' in path.parts:raise ValueError('unsafe ZIP entry')
                h=hashlib.sha256();size=0
                with z.open(name) as f:
                    for b in iter(lambda:f.read(1<<20),b''):h.update(b);size+=len(b)
                if h.hexdigest()!=expected[name]['sha256'] or size!=expected[name]['bytes']:raise ValueError('sealed entry differs: '+name)
            for path in sorted(HERE.glob('c01_regen_*.py')):
                entry=expected.get('source/'+path.name)
                if entry is None or digest(path)!=entry['sha256']:raise ValueError('diagnostic source differs: '+path.name)
    return dict(archive_verified=True,c01_passed=False,heavy_replay_executed=False)

if __name__=='__main__':print(json.dumps(verify()))
