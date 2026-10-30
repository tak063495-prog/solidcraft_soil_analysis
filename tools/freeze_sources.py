"""Freeze build, examples and verification sources before an acceptance run."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile


def manifest():
    root=Path(__file__).resolve().parents[1]
    files=[root/n for n in ('Cargo.toml','Cargo.lock','rust-toolchain.toml','requirements.lock.txt','sources.lock.json','README.md')]
    if (root/'requirements-reference.lock.txt').exists():files.append(root/'requirements-reference.lock.txt')
    if (root/'requirements-verification-win-cp313.lock.txt').exists():files.append(root/'requirements-verification-win-cp313.lock.txt')
    if (root/'.github').exists():
        files += [p for p in (root/'.github').rglob('*') if p.is_file() and p.suffix in ('.yml','.yaml')]
    files += [p for name in ('src','tests','examples','scripts','docs','tools') for p in (root/name).rglob('*')
              if p.is_file() and p.suffix in ('.rs','.py','.json','.ps1','.md','.csv')]
    return {p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(files))}


def freeze(output):
    root=Path(__file__).resolve().parents[1];sources=manifest()
    assert not (output/'source-snapshot.json').exists() and not (output/'source-snapshot.zip').exists()
    with zipfile.ZipFile(output/'source-snapshot.zip','x',compression=zipfile.ZIP_DEFLATED) as archive:
        for name,digest in sources.items():
            data=(root/name).read_bytes();assert hashlib.sha256(data).hexdigest()==digest;archive.writestr(name,data)
    with (output/'source-snapshot.json').open('x',encoding='utf-8') as stream:json.dump(sources,stream,indent=2)
    print('Frozen acceptance sources:',len(sources))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();freeze(a.output.resolve())
