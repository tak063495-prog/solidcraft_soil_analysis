"""Stream verified evidence into bounded archive parts, without a joined copy."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


class PartsWriter:
    def __init__(self, folder):
        self.folder = folder; self.parts = []; self.stream = None
        self.total = hashlib.sha256(); self.total_size = 0

    def write(self, raw):
        original_size = len(raw)
        self.total.update(raw); self.total_size += original_size
        while raw:
            if self.stream is None:
                self.name = f'evidence.tar.gz.part{len(self.parts):03d}'
                self.stream = (self.folder/self.name).open('xb')
                self.checksum = hashlib.sha256(); self.size = 0
            chunk, raw = raw[:16*1024**2-self.size], raw[16*1024**2-self.size:]
            self.stream.write(chunk); self.checksum.update(chunk); self.size += len(chunk)
            if self.size == 16*1024**2:
                self.finish_part()
        return original_size

    def finish_part(self):
        if self.stream is not None:
            self.stream.close(); self.stream = None
            self.parts.append(dict(name=self.name, size=self.size, sha256=self.checksum.hexdigest()))

    def flush(self):
        if self.stream is not None:
            self.stream.flush()


def pack(evidence, output):
    files = {}; omitted_empty_logs = 0
    for path in evidence.rglob('*'):
        if path.is_file() and 'observer-target' not in path.relative_to(evidence).parts and '__pycache__' not in path.parts:
            if path.suffix == '.log' and path.stat().st_size == 0:
                omitted_empty_logs += 1
                continue
            files[path.relative_to(evidence).as_posix()] = path
    source_paths = set(subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines())
    source_paths.update(p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py'))
    for name in sorted(source_paths):
        path = ROOT/name
        if path.suffix in ['.py','.rs','.toml','.lock','.yml'] or (name.startswith('examples/') and path.suffix=='.json'):
            files['sources/'+name] = path
    if any((output/p).exists() for p in ['archive-sha256.json','evidence.tar.gz.part000']):
        raise ValueError('archive output already exists')
    writer = PartsWriter(output); entries = {}
    with tarfile.open(fileobj=writer, mode='w|gz', compresslevel=6) as archive:
        for name, path in sorted(files.items()):
            if path.is_symlink():
                raise ValueError('symlink not allowed in evidence')
            info = tarfile.TarInfo(name); info.size = path.stat().st_size; info.mtime = 0; info.mode = 0o644
            with path.open('rb') as stream:
                checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
            with path.open('rb') as stream:
                archive.addfile(info, stream)
            entries[name] = dict(size=info.size, sha256=checksum)
    writer.finish_part()
    manifest = dict(parts=writer.parts, archive_sha256=writer.total.hexdigest(), archive_size=writer.total_size,
        entries=entries, all_material_contact_fields_preserved=True, all_new_linear_matrices_archived=False,
        omitted_empty_logs=omitted_empty_logs,
        new_linear_archive_scope='all SHA/rhs/solution/residual summaries and last 128 actual matrices per process')
    (output/'archive-sha256.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return dict(parts=len(writer.parts), entries=len(entries), archive_size=writer.total_size)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--output', type=Path, default=HERE); args=p.parse_args()
    print(json.dumps(pack(args.evidence.resolve(),args.output.resolve()),indent=2))
