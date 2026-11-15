"""Seal/verify new diagnostic evidence; verification never reruns equilibrium."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile
from a03follow_run import write


def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def verify(folder, extract=None):
    m=json.loads((folder/'archive-index.json').read_text(encoding='utf-8'))
    h=hashlib.sha256();size=0
    if extract:extract.mkdir(parents=True,exist_ok=False)
    with tempfile.TemporaryFile() as stream:
        for p in m['parts']:
            path=folder/p['name'];assert path.stat().st_size==p['size'] and sha(path)==p['sha256']
            with path.open('rb') as f:
                for raw in iter(lambda:f.read(1024*1024),b''):stream.write(raw);h.update(raw);size+=len(raw)
        assert size==m['archive_bytes'] and h.hexdigest()==m['archive_sha256'];stream.seek(0)
        with zipfile.ZipFile(stream) as z:
            assert len(z.namelist())==len(set(z.namelist())) and set(z.namelist())==set(m['entries'])
            for name in z.namelist():
                path=Path(name);assert not path.is_absolute() and '..' not in path.parts
                raw=z.read(name);e=m['entries'][name]
                assert len(raw)==e['size'] and hashlib.sha256(raw).hexdigest()==e['sha256'],name
                if extract:
                    target=extract/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
    return dict(parts=len(m['parts']),entries=len(m['entries']),all_sha256_verified=True,integrity_only_not_fresh_solve=True)


def seal(run, destination, retained):
    entries={};archive=run.parent/'a03follow-new-evidence.zip'
    sources=[('formal',run)]+[('retained/'+p.name,p) for p in retained]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for prefix,folder in sources:
            for path in sorted(folder.rglob('*')):
                if not path.is_file() or '__pycache__' in path.parts:continue
                name=prefix+'/'+path.relative_to(folder).as_posix();assert name not in entries
                z.write(path,name);entries[name]=dict(size=path.stat().st_size,sha256=sha(path))
    parts=[]
    with archive.open('rb') as f:
        while raw:=f.read(24*1024*1024):
            path=destination/('evidence.zip.part'+str(len(parts)).zfill(3))
            with path.open('xb') as part:part.write(raw)
            parts.append(dict(name=path.name,size=len(raw),sha256=sha(path)))
    write(destination/'archive-index.json',dict(archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),
        parts=parts,entries=entries,raw_inputs_states_checkpoints_refusals_and_failures_retained=True,
        solver_binaries_archived=False))
    return verify(destination)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',type=Path);p.add_argument('--extract',type=Path)
    p.add_argument('--run',type=Path);p.add_argument('--destination',type=Path);p.add_argument('--retained',type=Path,action='append',default=[])
    a=p.parse_args()
    print(verify(a.verify,a.extract) if a.verify else seal(a.run,a.destination,a.retained))


if __name__=='__main__':main()
