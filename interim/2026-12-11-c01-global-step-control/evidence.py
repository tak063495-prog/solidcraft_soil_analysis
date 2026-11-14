"""Strict archive integrity and field remeasurement; never launches a numerical solve."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import zipfile
HERE=Path(__file__).resolve().parent

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def extract(output):
    acceptance=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
    if acceptance.get('passed') is not False or acceptance.get('archive_index_sha256')!=digest(HERE/'archive-index.json'): raise ValueError('pinned archive index integrity')
    index=json.loads((HERE/'archive-index.json').read_text(encoding='utf-8'))
    output.mkdir(parents=True,exist_ok=False)
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        archive=Path(tmp)/'evidence.zip'
        with archive.open('xb') as target:
            for part in index['parts']:
                path=HERE/part['file']
                if path.stat().st_size!=part['bytes'] or digest(path)!=part['sha256']: raise ValueError('archive part integrity')
                with path.open('rb') as f: shutil.copyfileobj(f,target)
        if digest(archive)!=index['sha256']: raise ValueError('whole archive integrity')
        expected={e['path']:e for e in index['entries']}
        with zipfile.ZipFile(archive) as z:
            if len(z.namelist())!=len(expected) or set(z.namelist())!=set(expected): raise ValueError('archive inventory')
            for name in z.namelist():
                path=Path(name)
                if path.is_absolute() or '..' in path.parts: raise ValueError('unsafe archive path')
                target=output/path;target.parent.mkdir(parents=True,exist_ok=True)
                with z.open(name) as f,target.open('xb') as to: shutil.copyfileobj(f,to)
                if target.stat().st_size!=expected[name]['bytes'] or digest(target)!=expected[name]['sha256']: raise ValueError('entry integrity: '+name)
    return output

def replay(folder):
    spec=importlib.util.spec_from_file_location('c01_global_measure_replay',HERE/'measure.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    report,events=module.measure(folder/'paths')
    expected=json.loads((HERE/'measurements.json').read_text(encoding='utf-8'))
    module.validate_report(expected)
    if report!=expected or events!=json.loads((HERE/'switch-events.json').read_text(encoding='utf-8')): raise ValueError('remeasured fields differ')
    return dict(archive_verified=True,measurements_exact=True,numerical_acceptance=False)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    print(json.dumps(replay(extract(a.output)),ensure_ascii=False))
if __name__=='__main__': main()
