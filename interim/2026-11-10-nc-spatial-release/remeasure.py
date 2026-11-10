"""Verify all archive bytes and regenerate measurements without solving MCC."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import importlib.util
import gzip
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location('spatial_evidence_measure', HERE/'measure.py')
measure = importlib.util.module_from_spec(spec)
spec.loader.exec_module(measure)


def isolated_command(*arguments):
    return [sys.executable, '-I', '-X', 'utf8', *map(str, arguments)]


def complete_snapshot(snapshot):
    """Fill omitted imports only with bytes matching the original fingerprints.

    The archive already seals its diagnostic modules. Its original fingerprint
    file pins remaining dependencies. Never import a newer live JOINT module.
    """
    fingerprints = json.loads((HERE/'source-fingerprints.json').read_text())
    retained = ROOT/'interim/2026-11-11-joint-crack-tip/historical-spatial-dependencies'
    copied = {}
    for name, expected in fingerprints.items():
        path = Path(name)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('unsafe dependency path')
        if path.suffix not in ('.py', '.json', '.csv') or (snapshot/path).exists():
            continue
        candidates = [ROOT/path]
        if path.parent == Path('tools'):
            candidates.append(retained/path.name)
        for candidate in candidates:
            if not candidate.is_file():
                continue
            raw = candidate.read_bytes()
            variants = [raw, raw.replace(b'\r\n', b'\n'), raw.replace(b'\r\n', b'\n').replace(b'\n', b'\r\n')]
            payload = next((b for b in variants if hashlib.sha256(b).hexdigest() == expected), None)
            if payload is not None:
                target = snapshot/path; target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(payload)
                copied[name] = dict(sha256=expected, source=str(candidate.relative_to(ROOT)))
                break
    return copied


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
    evidence = out/'evidence'; seen = set(); eol_differences=[]; checkout_differences=[]
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
                current=ROOT/Path(*path.parts[1:])
                checkout=current.read_bytes() if current.is_file() else None
                if raw!=checkout:
                    difference=dict(path=ent.name,snapshot_sha256=hashlib.sha256(raw).hexdigest(),checkout_sha256=hashlib.sha256(checkout).hexdigest() if checkout is not None else None)
                    if checkout is not None and raw.replace(b'\r\n',b'\n')==checkout.replace(b'\r\n',b'\n'):
                        eol_differences.append(dict(difference,only_CRLF_LF_differs=True))
                    else:
                        checkout_differences.append(difference)
    reader.close()
    if seen != set(manifest['entries']):
        raise ValueError('archive inventory differs')
    snapshot=evidence/'sources'
    dependencies=complete_snapshot(snapshot)
    script=snapshot/'interim/2026-11-10-nc-spatial-release/measure.py'
    if not script.is_file():raise ValueError('sealed measurement source missing')
    child_code="""import importlib.util,json,sys
from pathlib import Path
script,evidence,output=map(Path,sys.argv[1:4])
spec=importlib.util.spec_from_file_location('sealed_spatial_measure',script)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
output.write_text(json.dumps(module.measure(evidence),indent=2,allow_nan=False)+'\\n',encoding='utf-8')
"""
    subprocess.run(isolated_command('-c',child_code,script,evidence,out/'measurements.json'),cwd=snapshot,check=True)
    actual=(out/'measurements.json').read_bytes()
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
    result = dict(archive_verified=True, verified_entries=len(seen), measurements_byte_identical=identical,discrete_evidence_exact=True,derived_float_differences=differences,numerical_acceptance=False,checkout_source_EOL_differences=eol_differences,checkout_source_differences=checkout_differences,pinned_missing_dependencies=dependencies,measurement_source='SHA-verified sealed archive snapshot in isolated subprocess',current_checkout_numerical_certification=False,jsonl_extracted_as_gzip=True,joined_archive_created=False)
    (out/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--allow-platform-float-differences',action='store_true')
    args=p.parse_args()
    print(json.dumps(run(args.output.resolve(),args.allow_platform_float_differences), indent=2))
