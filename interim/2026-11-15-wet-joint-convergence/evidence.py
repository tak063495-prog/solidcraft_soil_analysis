"""SHA-bound archive extraction and exact remeasurement; no new solve or certificate."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import zipfile
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from study import read, write, sha


def extract(destination):
    manifest=read(HERE/'archive-sha256.json')
    destination.mkdir(parents=True,exist_ok=False)
    archive=destination/'evidence.zip'
    with archive.open('xb') as f:
        for part in manifest['parts']:
            path=HERE/part['name']
            assert path.stat().st_size==part['bytes'] and sha(path)==part['sha256'], 'changed archive part'
            with path.open('rb') as src:shutil.copyfileobj(src,f)
    assert sha(archive)==manifest['archive_sha256']
    with zipfile.ZipFile(archive) as z:
        assert len(z.namelist())==len(set(z.namelist()))
        assert set(z.namelist())==set(manifest['entries']), 'changed entry inventory'
        for name,digest in manifest['entries'].items():
            p=Path(name)
            assert not p.is_absolute() and '..' not in p.parts
            data=z.read(name); assert hashlib.sha256(data).hexdigest()==digest, ('changed entry',name)
            path=destination/p;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    archive.unlink()
    return manifest


def validate_report(value):
    for flag in ('passed','all_history_accepted','p6_complete'):
        if value[flag] is not False:
            raise ValueError('unsupported all-history acceptance flag: '+flag)
    for row in value['space_only'].values():
        if row['spatial_convergence_accepted'] is not False or row['pointwise_contact_fields_compared'] is not False:
            raise ValueError('face-mean diagnostic falsely certifies pointwise convergence')
    for row in value['event_aligned_diagnostic'].values():
        if row['accepted'] is not False or row['existing_metrics_replaced'] is not False:
            raise ValueError('aligned diagnostic replaces existing metric')
    return value


def replay(destination):
    manifest=extract(destination)
    from measure import measure
    expected=validate_report(read(destination/'measurements.json'))
    actual,events=measure(destination/'paths');validate_report(actual)
    if actual!=expected or events!=read(destination/'switch-events.json'):
        raise ValueError('remeasured fields, events or acceptance flags differ')
    return dict(all_entries_sha256_verified=len(manifest['entries']),all_measurements_exact=True,
                software_replay_passed=True,numerical_acceptance=False,solver_rerun=False,
                source='Current diagnostic code reading SHA-bound immutable saved fields; sealed source files retained for provenance.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=replay(a.output.resolve());write(a.output/'replay-verification.json',result)
    print(json.dumps(result))
