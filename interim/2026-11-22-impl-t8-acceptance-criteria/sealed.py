"""Read a selected member without unpacking other meshes or modifying the seal."""
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import sys
import zipfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-12-17-c02-spatial-cause-diagnosis"))
from c02diag_extract import SplitZip


@contextmanager
def archive(folder, index_name):
    index = json.loads((folder / index_name).read_bytes())
    total = hashlib.sha256()
    parts = []
    for part in index["parts"]:
        path = folder / (part.get("name") or part["file"])
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            while block := stream.read(1024 * 1024):
                digest.update(block)
                total.update(block)
        if digest.hexdigest() != part["sha256"]:
            raise ValueError(f"invalid part {path}")
        parts.append(path)
    expected = index.get("archive_sha256") or index.get("zip_sha256") or index["sha256"]
    if total.hexdigest() != expected:
        raise ValueError(f"invalid combined archive {folder}")
    with SplitZip(parts) as raw, zipfile.ZipFile(raw) as z:
        yield z
