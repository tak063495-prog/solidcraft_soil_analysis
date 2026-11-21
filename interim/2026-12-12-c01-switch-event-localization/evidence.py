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
    compatibility=json.loads((HERE/'compatibility.json').read_text(encoding='utf-8'))
    if compatibility!=json.loads((folder/'compatibility/compatibility.json').read_text(encoding='utf-8')): raise ValueError('compatibility record differs')
    if digest(HERE/'baseline-source-tree.json')!=compatibility['baseline_blob_manifest_sha256']: raise ValueError('baseline source binding changed')
    prior_index=HERE.parents[0]/'2026-12-11-c01-global-step-control/archive-index.json'
    if digest(prior_index)!=compatibility['prior_index_sha256']: raise ValueError('prior evidence index changed')
    prior_entries={e['path']:e for e in json.loads(prior_index.read_text(encoding='utf-8'))['entries']}
    spec=importlib.util.spec_from_file_location('c01_path_compat_replay',HERE/'compatibility.py')
    compat=importlib.util.module_from_spec(spec);spec.loader.exec_module(compat)
    for row in compatibility['rows']:
        raw=folder/'compatibility'/row['mode']/f"g{row['steps']}"/'raw'
        if digest(raw/'accepted.jsonl')!=row['accepted_sha256'] or row['accepted_sha256']!=row['expected_sha256']: raise ValueError('legacy committed bytes changed')
        if row['mode']=='switch':
            expected=prior_entries[f"paths/g{row['steps']}/native-raw/accepted.jsonl"]['sha256']
            old=folder/'compatibility/prior-run'/f"g{row['steps']}"/'run.json'
            if digest(old)!=prior_entries[f"paths/g{row['steps']}/native-raw/run.json"]['sha256'] or expected!=row['expected_sha256']: raise ValueError('prior switch source changed')
        else:
            old=folder/'compatibility/default-base/run.json'
            if digest(folder/'compatibility/default-base/accepted.jsonl')!=row['expected_sha256']: raise ValueError('default baseline changed')
        if compat.canonical_bytes(json.loads(old.read_text(encoding='utf-8')))!=compat.canonical_bytes(json.loads((raw/'run.json').read_text(encoding='utf-8'))): raise ValueError('legacy full Run changed')
    retained=json.loads((HERE/'retained-failures/index.json').read_text(encoding='utf-8'))
    archive=HERE/'retained-failures/prototypes.zip'
    if digest(archive)!=retained['archive_sha256']: raise ValueError('prototype failure archive changed')
    with zipfile.ZipFile(archive) as z:
        if set(z.namelist())!={e['path'] for e in retained['entries']}: raise ValueError('prototype inventory changed')
        for e in retained['entries']:
            data=z.read(e['path'])
            if len(data)!=e['bytes'] or hashlib.sha256(data).hexdigest()!=e['sha256']: raise ValueError('prototype bytes changed')
    return dict(archive_verified=True,measurements_exact=True,legacy_accepted_bytes_exact=True,
                legacy_runs_exact_except_clock_and_executable=True,retained_failures_verified=True,numerical_acceptance=False)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    print(json.dumps(replay(extract(a.output)),ensure_ascii=False))
if __name__=='__main__': main()
