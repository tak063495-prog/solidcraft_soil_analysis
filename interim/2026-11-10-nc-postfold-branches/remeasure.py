"""Verify sealed bytes and independently reassemble fields; no new equilibrium solve."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile

HERE=Path(__file__).resolve().parent


class Reader(io.RawIOBase):
    def __init__(self,paths):self.paths=iter(paths);self.stream=None
    def readable(self):return True
    def read(self,size=-1):
        if size<0:raise ValueError('bounded stream reads required')
        data=bytearray()
        while len(data)<size:
            if self.stream is None:
                try:self.stream=next(self.paths).open('rb')
                except StopIteration:break
            block=self.stream.read(size-len(data))
            if block:data.extend(block)
            else:self.stream.close();self.stream=None
        return bytes(data)
    def close(self):
        if self.stream is not None:self.stream.close()
        super().close()


def run(out):
    out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((HERE/'archive-sha256.json').read_text(encoding='utf-8'))
    full=hashlib.sha256();total=0
    for part in manifest['parts']:
        sha=hashlib.sha256();size=0
        with (HERE/part['name']).open('rb') as stream:
            for block in iter(lambda:stream.read(1024*1024),b''):
                sha.update(block);full.update(block);size+=len(block);total+=len(block)
        if size!=part['size'] or sha.hexdigest()!=part['sha256']:raise ValueError('archive part differs')
    if total!=manifest['archive_size'] or full.hexdigest()!=manifest['archive_sha256']:raise ValueError('complete archive differs')
    seen=set();destination=out/'evidence'
    reader=Reader([HERE/p['name'] for p in manifest['parts']])
    try:
        with tarfile.open(fileobj=reader,mode='r|gz') as archive:
            for entry in archive:
                path=Path(entry.name)
                if not entry.isfile() or path.is_absolute() or '..' in path.parts or entry.name in seen or entry.name not in manifest['entries']:
                    raise ValueError('unsafe or unlisted archive entry')
                expected=manifest['entries'][entry.name]
                if entry.size!=expected['size']:raise ValueError('entry size differs')
                target=destination/path;target.parent.mkdir(parents=True,exist_ok=True)
                sha=hashlib.sha256()
                with archive.extractfile(entry) as source,target.open('xb') as stream:
                    for block in iter(lambda:source.read(1024*1024),b''):sha.update(block);stream.write(block)
                if sha.hexdigest()!=expected['sha256']:raise ValueError('entry sha differs')
                seen.add(entry.name)
    finally:reader.close()
    if seen!=set(manifest['entries']):raise ValueError('entry inventory differs')
    driver=destination/'drivers/measure.py'
    result=out/'measurements.json'
    code='import runpy,sys;sys.path.insert(0,'+repr(str(driver.parent))+');sys.argv=['+repr(str(driver))+',"--source",'+repr(str(destination/'study'))+',"--output",'+repr(str(result))+'];runpy.run_path('+repr(str(driver))+',run_name="__main__")'
    subprocess.run([sys.executable,'-I','-X','utf8','-c',code],check=True)
    if result.read_bytes()!=(HERE/'measurements.json').read_bytes():raise ValueError('sealed raw-evidence measurements differ')
    record=dict(archive_verified=True,verified_entries=len(seen),measurements_byte_identical=True,numerical_acceptance=False,
                prior_input_provenance='prepared/prior/extraction.json; independently verified prior archive at calculation time')
    (out/'verification.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
