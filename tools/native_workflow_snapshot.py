"""Freeze the external native analysis scripts without selecting an executable.

Runtime paths and interpreter/solver programs remain trusted host arguments.
Imported CAD files never supply script or executable paths.
"""
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import time
import uuid


def atomic_json(path, value):
    temporary = path.with_name(path.name+'.'+uuid.uuid4().hex+'.part')
    with temporary.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, allow_nan=False)
    # A Windows observer may briefly hold the old file without delete sharing.
    # Retry replacement of this same owned, complete temporary; never rewrite
    # the live file or select any other temporary/path.
    for attempt in range(5):
        try:
            os.replace(temporary, path)
            return
        except PermissionError:
            if attempt == 4:
                raise
            time.sleep(.01)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def capture(root, output, bindings):
    root, output = root.resolve(), output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for name, digest in bindings.items():
        relative = PurePosixPath(name)
        if relative.is_absolute() or '..' in relative.parts or '\\' in name or ':' in name:
            raise ValueError('workflow inventory must contain safe relative paths')
        source = (root/name).resolve()
        target = (output/name).resolve()
        if not source.is_relative_to(root) or not target.is_relative_to(output):
            raise ValueError('workflow source or destination escaped its root')
        if sha(source) != digest:
            raise ValueError('workflow source changed before snapshot: '+name)
        target.parent.mkdir(parents=True, exist_ok=True)
        with source.open('rb') as incoming, target.open('xb') as outgoing:
            for block in iter(lambda: incoming.read(1024*1024), b''):
                outgoing.write(block)
        if sha(source) != digest or sha(target) != digest:
            raise ValueError('workflow source changed during snapshot: '+name)
    verify(output, bindings)
    return output


def verify(root, bindings):
    root = root.resolve()
    for name, digest in bindings.items():
        path = (root/name).resolve()
        if not path.is_relative_to(root) or sha(path) != digest:
            raise ValueError('captured workflow source changed: '+name)
