"""Verify every new evidence byte; remeasure with isolated sealed sources."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent


class PartsReader:
    def __init__(self, paths):
        self.paths = iter(paths)
        self.current = None

    def read(self, size):
        data = bytearray()
        while len(data) < size:
            if self.current is None:
                try:
                    self.current = next(self.paths).open('rb')
                except StopIteration:
                    break
            block = self.current.read(size-len(data))
            if block:
                data.extend(block)
            else:
                self.current.close()
                self.current = None
        return bytes(data)

    def close(self):
        if self.current is not None:
            self.current.close()


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((HERE/'archive-sha256.json').read_text())
    total = hashlib.sha256()
    size = 0
    for part in manifest['parts']:
        h = hashlib.sha256()
        count = 0
        with (HERE/part['name']).open('rb') as stream:
            for block in iter(lambda: stream.read(1024*1024), b''):
                h.update(block)
                total.update(block)
                count += len(block)
                size += len(block)
        if h.hexdigest() != part['sha256'] or count != part['size']:
            raise ValueError('new evidence part changed')
    if total.hexdigest() != manifest['archive_sha256'] or size != manifest['archive_size']:
        raise ValueError('new complete archive changed')
    evidence = out/'evidence'
    seen = set()
    reader = PartsReader([HERE/p['name'] for p in manifest['parts']])
    try:
        with tarfile.open(fileobj=reader, mode='r|gz') as archive:
            for entry in archive:
                p = Path(entry.name)
                if not entry.isfile() or p.is_absolute() or '..' in p.parts or entry.name in seen or entry.name not in manifest['entries']:
                    raise ValueError('unsafe or unlisted new evidence member')
                expected = manifest['entries'][entry.name]
                if entry.size != expected['size']:
                    raise ValueError('new evidence member size differs')
                target = evidence/p
                target.parent.mkdir(parents=True, exist_ok=True)
                h = hashlib.sha256()
                with archive.extractfile(entry) as source, target.open('xb') as destination:
                    for block in iter(lambda: source.read(1024*1024), b''):
                        h.update(block)
                        destination.write(block)
                if h.hexdigest() != expected['sha256']:
                    raise ValueError('new evidence member changed')
                seen.add(entry.name)
    finally:
        reader.close()
    if seen != set(manifest['entries']):
        raise ValueError('new evidence inventory differs')
    snapshot = evidence/'sources'
    sealed_extract = snapshot/'interim/2026-11-10-nc-tangent-branch-diagnosis/extract.py'
    if sealed_extract.read_bytes() != (HERE/'extract.py').read_bytes():
        raise ValueError('archive verifier source differs from sealed copy')
    spec = importlib.util.spec_from_file_location('verified_prior_archive_extractor', HERE/'extract.py')
    extractor = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(extractor)
    prior = out/'previous-evidence'
    proof = extractor.run(prior)
    if proof != json.loads((evidence/'previous-extraction.json').read_text()):
        raise ValueError('previous evidence proof differs from sealed diagnosis input')
    script = snapshot/'interim/2026-11-10-nc-tangent-branch-diagnosis/measure.py'
    subprocess.run([sys.executable, '-I', '-X', 'utf8', str(script), '--source', str(evidence/'diagnosis'),
                    '--evidence', str(prior), '--output', str(out/'measurements.json')], check=True, cwd=snapshot)
    identical = (out/'measurements.json').read_bytes() == (HERE/'measurements.json').read_bytes()
    joint = json.loads((evidence/'joint-history/measurements.json').read_text())
    for count in (4, 8):
        folder = evidence/'joint-history'/f'g{count}'
        output = out/f'joint-g{count}.json'
        subprocess.run([sys.executable, '-I', '-X', 'utf8', str(snapshot/'tools/audit_wet_joint_accepted_history.py'),
                        '--history', str(folder/'accepted'), '--result', str(folder/'result.json'),
                        '--output', str(output)], check=True, cwd=snapshot)
        if json.loads(output.read_text()) != joint['checks'][str(count)]['audit']:
            raise ValueError('joint propagation remeasurement differs')
    result = dict(archive_verified=True, verified_entries=len(seen), previous_archive_verified=True,
                  measurements_byte_identical=identical, joint_propagation_checks_exact=True,
                  numerical_acceptance=False, full_release_accepted=False,
                  measurement_source='SHA-verified isolated sealed sources')
    (out/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    if not identical:
        raise ValueError('measurement JSON differs; record platform differences, never infer acceptance')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    print(json.dumps(run(p.parse_args().output.resolve()), indent=2))
