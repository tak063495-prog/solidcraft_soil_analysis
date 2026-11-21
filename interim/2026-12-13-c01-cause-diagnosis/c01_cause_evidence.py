"""Verify sealed diagnostic bytes and recompute measurements without any new solve."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile
import c01_cause_diagnostics as diag

HERE = Path(__file__).resolve().parent


def extract(output):
    acceptance = diag.read(HERE/'acceptance.json')
    if (acceptance['passed'] is not False or acceptance['archive_index_sha256'] != diag.sha(HERE/'archive-index.json')
            or acceptance['geometry_sha256'] != diag.sha(HERE/'geometry.json')):
        raise ValueError('diagnostic index/acceptance changed')
    index = diag.read(HERE/'archive-index.json')
    output.mkdir(parents=True,exist_ok=False)
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        archive = Path(tmp)/'evidence.zip'
        with archive.open('xb') as out:
            for part in index['parts']:
                path = HERE/part['file']
                if path.stat().st_size != part['bytes'] or diag.sha(path) != part['sha256']:
                    raise ValueError('diagnostic part integrity')
                with path.open('rb') as inp:
                    shutil.copyfileobj(inp,out)
        if diag.sha(archive) != index['sha256']:
            raise ValueError('diagnostic ZIP integrity')
        expected = {e['path']:e for e in index['entries']}
        with zipfile.ZipFile(archive) as z:
            if len(z.namelist()) != len(expected) or set(z.namelist()) != set(expected):
                raise ValueError('diagnostic inventory')
            for name in z.namelist():
                path = Path(name)
                if path.is_absolute() or '..' in path.parts:
                    raise ValueError('unsafe archive path')
                dest = output/path
                dest.parent.mkdir(parents=True,exist_ok=True)
                with z.open(name) as inp,dest.open('xb') as out:
                    shutil.copyfileobj(inp,out)
                if dest.stat().st_size != expected[name]['bytes'] or diag.sha(dest) != expected[name]['sha256']:
                    raise ValueError('diagnostic entry integrity: '+name)
    return output


def replay(output):
    data = extract(output/'diagnostic')
    for source in sorted(HERE.glob('c01_cause_*.py')):
        if diag.sha(source) != diag.sha(data/'final-source'/source.name):
            raise ValueError('diagnostic source bytes differ: '+source.name)
    old,_ = diag.dependencies()
    for control,folder in (('switch','2026-12-11-c01-global-step-control'),('path','2026-12-12-c01-switch-event-localization')):
        selected = output/('saved-'+control)
        diag.extract_selected(diag.ROOT/'interim'/folder,selected)
        actual,_ = diag.postprocess(selected/'paths',old,data/'native-counts'/control)
        expected = diag.read(HERE/f'measurement-{control}.json')
        if actual['pairs'] != expected['pairs'] or actual['audits'] != expected['audits']:
            raise ValueError('independent postprocessing differs: '+control)
    command = [sys.executable,str(HERE/'c01_cause_compare.py'),'--saved-switch',str(output/'saved-switch/paths'),'--saved-path',str(output/'saved-path/paths'),'--experiments',str(data/'experiments'),'--output',str(output/'comparison')]
    subprocess.run(command,check=True)
    if diag.read(output/'comparison/comparison.json') != diag.read(HERE/'comparison.json'):
        raise ValueError('reference/quadrature comparison differs')
    return dict(archive_verified=True,measurements_exact=True,numerical_solve=False,c01_accepted=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--integrity-only',action='store_true')
    args = p.parse_args()
    if args.integrity_only:
        extract(args.output)
        print(json.dumps(dict(archive_verified=True,numerical_solve=False)))
    else:
        args.output.mkdir(parents=True,exist_ok=False)
        print(json.dumps(replay(args.output)))


if __name__=='__main__':
    main()
