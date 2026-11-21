"""Authenticate sealed C-02 diagnostics without rerunning equilibrium or replay."""
import hashlib
import json
from pathlib import Path
import sys
import zipfile
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from c02diag_extract import SplitZip,digest


def document_status(path):
    # read_text performs universal newline decoding, including actual CRLF files.
    text=Path(path).read_text(encoding='utf-8')
    return text.rstrip('\n').endswith('C-02は未達、passed=false。')


def verify():
    index=json.loads((HERE/'archive-index.json').read_text(encoding='utf-8'))
    status=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
    if status['passed'] is not False or status['archive_index_sha256']!=digest(HERE/'archive-index.json'):
        raise ValueError('acceptance or archive index changed')
    parts=[];whole=hashlib.sha256()
    for item in index['parts']:
        p=HERE/item['file'];h=hashlib.sha256()
        if p.stat().st_size!=item['bytes']:raise ValueError('sealed part length differs')
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1<<20),b''):h.update(block);whole.update(block)
        if h.hexdigest()!=item['sha256']:raise ValueError('sealed part SHA differs')
        parts.append(p)
    if whole.hexdigest()!=index['sha256']:raise ValueError('sealed ZIP SHA differs')
    expected={x['path']:x for x in index['entries']}
    with SplitZip(parts) as raw,zipfile.ZipFile(raw) as archive:
        if len(archive.namelist())!=len(expected) or set(archive.namelist())!=set(expected):raise ValueError('sealed inventory differs')
        for name,item in sorted(expected.items()):
            p=Path(name)
            if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe ZIP entry')
            h=hashlib.sha256();size=0
            with archive.open(name) as f:
                for block in iter(lambda:f.read(1<<20),b''):h.update(block);size+=len(block)
            if h.hexdigest()!=item['sha256'] or size!=item['bytes']:raise ValueError('sealed entry differs: '+name)
        for p in sorted(HERE.glob('c02diag_*.py')):
            if expected.get('source/'+p.name,{}).get('sha256')!=digest(p):raise ValueError('diagnostic source differs: '+p.name)
    if status['new_equilibrium_paths_completed']!=8 or status['pc_executed'] is not False:
        raise ValueError('diagnostic scope changed')
    return dict(archive_verified=True,entries_verified=len(expected),c02_passed=False,heavy_recalculation_executed=False)

if __name__=='__main__':print(json.dumps(verify()))
