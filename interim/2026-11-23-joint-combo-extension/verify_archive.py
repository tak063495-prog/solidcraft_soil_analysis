"""Verify new split evidence bytes and every raw entry; no solving/rewriting."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile


def verify(folder):
    manifest=json.loads((folder/'evidence-manifest.json').read_text(encoding='utf-8'))
    joined=hashlib.sha256();size=0
    with tempfile.TemporaryFile() as archive:
        for part in manifest['parts']:
            path=folder/part['name'];h=hashlib.sha256();n=0
            with path.open('rb') as stream:
                for data in iter(lambda:stream.read(1024*1024),b''):
                    h.update(data);joined.update(data);archive.write(data);n+=len(data)
            assert n==part['size'] and h.hexdigest()==part['sha256'],('archive part differs',path.name)
            size+=n
        assert size==manifest['archive_bytes'] and joined.hexdigest()==manifest['archive_sha256'],'joined archive differs'
        archive.seek(0)
        with zipfile.ZipFile(archive) as z:
            names=z.namelist();assert len(names)==len(set(names)) and set(names)==set(manifest['entries']),'inventory differs'
            for name in names:
                path=Path(name);assert not path.is_absolute() and '..' not in path.parts
                h=hashlib.sha256();n=0
                with z.open(name) as stream:
                    for data in iter(lambda:stream.read(1024*1024),b''):h.update(data);n+=len(data)
                expected=manifest['entries'][name]
                assert n==expected['size'] and h.hexdigest()==expected['sha256'],('entry differs',name)
    return dict(archive_bytes_verified=True,every_entry_verified=True,parts=len(manifest['parts']),entries=len(names),
                archive_sha256=joined.hexdigest(),saved_evidence_audit_only=True,counted_as_fresh_solve=False)


def main():
    p=argparse.ArgumentParser();p.add_argument('folder',type=Path);p.add_argument('--report',type=Path);a=p.parse_args()
    result=verify(a.folder)
    if a.report:
        with a.report.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
