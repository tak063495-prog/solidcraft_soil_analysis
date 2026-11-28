"""Exclusive UTF-8/LF evidence writes and streaming fingerprints."""
import hashlib
import json
from pathlib import Path

def write(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()
