"""Losslessly compress diagnostic JSONL and split a checksummed evidence tar."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tarfile
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from freeze_sources import manifest as source_manifest
from linear_evidence_codec import encode, json_lines


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def compress(evidence):
    p = evidence/'compressed-jsonl.json'
    records = json.loads(p.read_text()) if p.exists() else {}
    for folder in sorted(evidence.iterdir()):
        if not folder.is_dir() or not ((folder/'report.json').exists() or (folder/'external-watchdog.json').exists()):
            continue
        for raw in sorted(folder.glob('*.jsonl')):
            q = Path(str(raw)+'.gz'); digest = sha(raw); size = raw.stat().st_size
            with raw.open('rb') as a, q.open('xb') as b, gzip.GzipFile(fileobj=b, filename='', mode='wb', compresslevel=1, mtime=0) as stream:
                for chunk in iter(lambda: a.read(1024*1024), b''):
                    stream.write(chunk)
            check = hashlib.sha256(); length = 0
            with gzip.open(q, 'rb') as stream:
                for chunk in iter(lambda: stream.read(1024*1024), b''):
                    check.update(chunk); length += len(chunk)
            if check.hexdigest() != digest or length != size:
                raise ValueError('compression round trip differs')
            records[raw.relative_to(evidence).as_posix()] = dict(raw_sha256=digest, raw_size=size, gzip_sha256=sha(q), gzip_size=q.stat().st_size)
            raw.unlink()
    for name,record in records.items():
        if not name.endswith('/linear-systems.jsonl') or record.get('encoding')=='float64-lossless':continue
        source=evidence/(name+'.gz');target=evidence/(name+'.f64.gz')
        encode(source,target);checksum=hashlib.sha256();length=0
        for line in json_lines(target):checksum.update(line);length+=len(line)
        if checksum.hexdigest()!=record['raw_sha256'] or length!=record['raw_size']:
            raise ValueError('float64 evidence did not reconstruct exact original JSON bytes')
        record.update(encoding='float64-lossless',encoded_sha256=sha(target),encoded_size=target.stat().st_size,bitwise_delta=True)
        source.unlink()
        p.write_text(json.dumps(records,indent=2)+'\n')
    p.write_text(json.dumps(records, indent=2)+'\n')


class SplitWriter:
    def __init__(self, output, size=16*1024*1024):
        self.output = output; self.size = size; self.parts = []; self.stream = None; self.position = 0
    def write(self, data):
        start = 0
        while start < len(data):
            if self.stream is None:
                p = self.output/f'evidence.tar.gz.part{len(self.parts):03d}'
                self.stream = p.open('xb'); self.parts.append(p); self.position = 0
            count = min(len(data)-start, self.size-self.position)
            self.stream.write(data[start:start+count]); self.position += count; start += count
            if self.position == self.size:
                self.stream.close(); self.stream = None
        return len(data)
    def flush(self):
        if self.stream:
            self.stream.flush()
    def close(self):
        if self.stream:
            self.stream.close(); self.stream = None


def package(evidence, output):
    compress(evidence)
    sources = source_manifest()
    capture=HERE.parent/'2026-11-09-nc-pseudo-arclength/capture_reference_history.py'
    sources[capture.relative_to(ROOT).as_posix()]=sha(capture)
    for path in HERE.glob('*.py'):
        sources[path.relative_to(ROOT).as_posix()] = sha(path)
    for p in evidence.glob('*/controls.json'):
        for name, expected in json.loads(p.read_text())['binding']['driver_sources'].items():
            if sha(ROOT/name) != expected:
                raise ValueError('numerical source changed since execution: '+name)
            sources[name] = expected
    for name, expected in sorted(sources.items()):
        q = evidence/'sources'/name; q.parent.mkdir(parents=True, exist_ok=True)
        q.write_bytes((ROOT/name).read_bytes())
        if sha(q) != expected:
            raise ValueError('source snapshot differs')
    (evidence/'source-snapshot.json').write_text(json.dumps(sources, indent=2)+'\n')
    files = []
    for folder in sorted(evidence.iterdir()):
        if not folder.is_dir():
            if folder.suffix in ('.json', '.log', '.txt'):
                files.append(folder)
            continue
        if folder.name in ('sources', 'fresh-native-prefix', 'fresh-prefix'):
            files += [p for p in folder.rglob('*') if p.is_file()]
        elif folder.name == 'prepared':
            files += [p for p in folder.rglob('*') if p.is_file() and 'observer-target' not in p.parts]
        elif (folder/'controls.json').exists() or (folder/'post-full-release.json').exists():
            files += [p for p in folder.rglob('*') if p.is_file() and (p.suffix != '.log' or p.stat().st_size)]
    files = sorted(set(files))
    entries = {}
    writer = SplitWriter(output)
    with tarfile.open(fileobj=writer, mode='w|gz', compresslevel=6) as archive:
        for p in files:
            name = p.relative_to(evidence).as_posix(); entries[name] = dict(sha256=sha(p), size=p.stat().st_size)
            archive.add(p, arcname=name, recursive=False)
    writer.close()
    total = hashlib.sha256(); parts = []
    for p in writer.parts:
        data = p.read_bytes(); total.update(data)
        parts.append(dict(name=p.name, size=len(data), sha256=hashlib.sha256(data).hexdigest()))
    result = dict(archive_sha256=total.hexdigest(), archive_size=sum(p['size'] for p in parts), parts=parts, entries=entries,
                  contents='All committed/full and tentative rate fields, Newton/rate trials, actual linear systems, checkpoints, own input histories, observer and diagnostic sources. Transient material request/response bytes replaced by their recorded hashes; executable binaries and empty material stdout logs omitted.')
    (output/'archive-sha256.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--compress-only', action='store_true')
    args = p.parse_args()
    if args.compress_only:
        compress(args.evidence.resolve())
    else:
        package(args.evidence.resolve(), args.output.resolve())
