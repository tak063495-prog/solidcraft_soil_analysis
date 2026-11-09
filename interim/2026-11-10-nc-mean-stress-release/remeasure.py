"""Verify all archive bytes and regenerate measurements without solving MCC."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import importlib.util
import gzip
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from linear_evidence_codec import json_lines

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location('mean_stress_evidence_measure', HERE/'measure.py')
measure = importlib.util.module_from_spec(spec)
spec.loader.exec_module(measure)


class PartsReader:
    """Sequential tar input without a second full-size joined archive."""
    def __init__(self, paths):
        self.paths=iter(paths); self.current=None
    def read(self, size):
        data=bytearray()
        while len(data)<size:
            if self.current is None:
                try:self.current=next(self.paths).open('rb')
                except StopIteration:break
            chunk=self.current.read(size-len(data))
            if chunk:data.extend(chunk)
            else:self.current.close();self.current=None
        return bytes(data)
    def close(self):
        if self.current is not None:self.current.close()


def run(out,allow_platform_float_differences=False):
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((HERE/'archive-sha256.json').read_text())
    total = hashlib.sha256(); size = 0
    for part in manifest['parts']:
        p = HERE/part['name']; raw = p.read_bytes()
        if len(raw) != part['size'] or hashlib.sha256(raw).hexdigest() != part['sha256']:
            raise ValueError('archive part changed: '+p.name)
        total.update(raw); size += len(raw)
    if total.hexdigest() != manifest['archive_sha256'] or size != manifest['archive_size']:
        raise ValueError('joined archive differs')
    evidence = out/'evidence'; seen = set(); eol_differences=[]
    reader=PartsReader([HERE/p['name'] for p in manifest['parts']])
    with tarfile.open(fileobj=reader,mode='r|gz') as z:
        for ent in z:
            path = Path(ent.name)
            if not ent.isfile() or path.is_absolute() or '..' in path.parts or ent.name in seen or ent.name not in manifest['entries']:
                raise ValueError('unsafe/unlisted archive entry')
            expected = manifest['entries'][ent.name]
            if ent.size!=expected['size']:raise ValueError('archive entry size differs')
            q=evidence/path;q.parent.mkdir(parents=True,exist_ok=True)
            compressed=path.suffix=='.jsonl'
            target=gzip.open(str(q)+'.gz','wb',compresslevel=1) if compressed else q.open('xb')
            checksum=hashlib.sha256();length=0
            with target as f:
                source=z.extractfile(ent)
                for chunk in iter(lambda:source.read(1024*1024),b''):
                    checksum.update(chunk);length+=len(chunk);f.write(chunk)
            if length != expected['size'] or checksum.hexdigest() != expected['sha256']:
                raise ValueError('archive entry differs: '+ent.name)
            seen.add(ent.name)
            if path.parts[0] == 'sources':
                raw=q.read_bytes()
                checkout=(ROOT/Path(*path.parts[1:])).read_bytes()
                if raw!=checkout:
                    if raw.replace(b'\r\n',b'\n')!=checkout.replace(b'\r\n',b'\n'):
                        raise ValueError('diagnostic source differs: '+ent.name)
                    eol_differences.append(dict(path=ent.name,snapshot_sha256=hashlib.sha256(raw).hexdigest(),checkout_sha256=hashlib.sha256(checkout).hexdigest(),only_CRLF_LF_differs=True))
    reader.close()
    if seen != set(manifest['entries']):
        raise ValueError('archive inventory differs')
    compressed=json.loads((evidence/'compressed-jsonl.json').read_text())
    for name,record in compressed.items():
        checksum=hashlib.sha256();length=0
        if record.get('encoding')=='float64-lossless':
            for chunk in json_lines(evidence/(name+'.f64.gz')):
                checksum.update(chunk);length+=len(chunk)
        else:
            with gzip.open(evidence/(name+'.gz'),'rb') as f:
                for chunk in iter(lambda:f.read(1024*1024),b''):
                    checksum.update(chunk);length+=len(chunk)
        if checksum.hexdigest()!=record['raw_sha256'] or length!=record['raw_size']:
            raise ValueError('compressed raw JSONL differs: '+name)
    actual = (json.dumps(measure.measure(evidence,allow_platform_float_differences), indent=2, allow_nan=False)+'\n').encode()
    (out/'measurements.json').write_bytes(actual)
    expected=(HERE/'measurements.json').read_bytes();identical=actual==expected;differences=[]
    def structure(a,b,path=''):
        if isinstance(a,float) and isinstance(b,float):
            if a!=b:differences.append(dict(path=path,recorded=a,remeasured=b,absolute_difference=abs(a-b)))
        elif type(a)!=type(b):raise ValueError('remeasured structure differs: '+path)
        elif isinstance(a,dict):
            if a.keys()!=b.keys():raise ValueError('remeasured keys differ: '+path)
            for key in a:structure(a[key],b[key],path+'/'+key)
        elif isinstance(a,list):
            if len(a)!=len(b):raise ValueError('remeasured count differs: '+path)
            for i,(aa,bb) in enumerate(zip(a,b)):structure(aa,bb,path+'/'+str(i))
        elif a!=b:raise ValueError('remeasured discrete evidence differs: '+path)
    structure(json.loads(expected),json.loads(actual))
    if not identical and not allow_platform_float_differences:
        raise ValueError('remeasured JSON differs')
    result = dict(archive_verified=True, verified_entries=len(seen), measurements_byte_identical=identical,discrete_evidence_exact=True,derived_float_differences=differences,numerical_acceptance=False,checkout_source_EOL_differences=eol_differences,jsonl_extracted_as_gzip=True,joined_archive_created=False)
    (out/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--allow-platform-float-differences',action='store_true')
    args=p.parse_args()
    print(json.dumps(run(args.output.resolve(),args.allow_platform_float_differences), indent=2))
