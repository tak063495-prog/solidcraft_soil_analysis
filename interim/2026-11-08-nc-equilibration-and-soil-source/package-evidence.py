"""Pack reproducible raw evidence and verify every archive member on extraction."""
import hashlib
import json
from pathlib import Path
import zipfile
out = Path(__file__).resolve().parent
root = out.parents[1]
folder = root/'work/evidence'
archive = out/'evidence.zip'
entries = {}
with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for path in sorted(folder.rglob('*')):
        if path.is_file():
            name = path.relative_to(folder).as_posix()
            entries[name] = hashlib.sha256(path.read_bytes()).hexdigest()
            method = zipfile.ZIP_LZMA if name in ('nc/dd/trace.jsonl', 'nc/accurate/trace.jsonl') else zipfile.ZIP_DEFLATED
            z.write(path, name, compress_type=method)
manifest = dict(archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(), entries_sha256=entries)
(out/'archive-sha256.json').write_text(json.dumps(manifest, indent=2)+'\n')
extracted = root/'work/reextracted-evidence-lzma'
extracted.mkdir(exist_ok=False)
with zipfile.ZipFile(archive) as z:
    assert set(z.namelist()) == set(entries)
    for name, digest in entries.items():
        data = z.read(name); assert hashlib.sha256(data).hexdigest() == digest
    z.extractall(extracted)
print('Verified archive entries:', len(entries), 'bytes:', archive.stat().st_size)
