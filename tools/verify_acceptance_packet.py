"""Read-only hash and source archive verification of a sealed acceptance run."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def verify(run):
    root=Path(__file__).resolve().parents[1]
    manifest=run/'final-validation.json'
    data=json.loads(manifest.read_text(encoding='utf-8-sig'))
    assert data['passed'] and data['schema'].startswith('soil-fem-')
    for name,expected in data['sha256'].items():
        path=(root/name).resolve()
        assert path.is_relative_to(root),'evidence path escapes workspace'
        assert sha(path)==expected,path
    archive=run/'source-snapshot.zip'
    assert sha(archive)==data['source_snapshot_sha256']
    sources=json.loads((run/'source-snapshot.json').read_text(encoding='utf-8-sig'))
    with zipfile.ZipFile(archive) as zipped:
        names={name.replace('\\','/'):name for name in zipped.namelist()}
        expected_sources={name.replace('\\','/'):value for name,value in sources.items()}
        assert len(zipped.namelist())==len(names) and len(sources)==len(expected_sources)
        assert set(names)==set(expected_sources)
        assert zipped.testzip() is None
        for name,expected in expected_sources.items():
            assert hashlib.sha256(zipped.read(names[name])).hexdigest()==expected,name
    assert sha(run/'solver/soil-fem-3d.exe')==data['solver_sha256']
    return {'schema':'soil-fem-acceptance-packet-integrity/1','passed':True,
        'manifest_sha256':sha(manifest),'evidence_files':len(data['sha256']),
        'frozen_source_files':len(sources),'solver_sha256':data['solver_sha256'],
        'scope':'Stored evidence and frozen sources; does not certify the current working tree or additional physics.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--output',type=Path)
    a=p.parse_args();report=verify(a.run.resolve())
    if a.output:
        with a.output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(report,indent=2))
