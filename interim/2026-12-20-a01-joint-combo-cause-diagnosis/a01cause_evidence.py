"""Verify/extract immutable evidence and seal only new observations, with SHA."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile
from a01cause_driver import write


def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def extract_old(folder,out):
    m=json.loads((folder/'evidence-manifest.json').read_text(encoding='utf-8'));out.mkdir(parents=True,exist_ok=False)
    digest=hashlib.sha256();count=0
    with tempfile.TemporaryFile() as stream:
        for part in m['parts']:
            path=folder/part['name'];assert path.stat().st_size==part['size'] and sha(path)==part['sha256']
            with path.open('rb') as f:
                for raw in iter(lambda:f.read(1024*1024),b''):stream.write(raw);digest.update(raw)
        assert digest.hexdigest()==m['archive_sha256'];stream.seek(0)
        with zipfile.ZipFile(stream) as z:
            assert len(z.namelist())==len(set(z.namelist())) and set(z.namelist())==set(m['entries'])
            for name in z.namelist():
                path=Path(name);assert not path.is_absolute() and '..' not in path.parts
                raw=z.read(name);entry=m['entries'][name]
                assert len(raw)==entry['size'] and hashlib.sha256(raw).hexdigest()==entry['sha256'];count+=1
                if name.startswith(('selected/','interacting-stacked-attempt/')) and 'restarts/' not in name and 'checkpoint-' not in name:
                    target=out/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
    result=dict(parts=len(m['parts']),entries=count,archive_sha256=digest.hexdigest(),all_entries_verified=True)
    write(out/'verified.json',result);return result


def seal(sources,destination):
    archive=sources[0][1].parent/('a01cause-new-'+destination.name+'.zip')
    entries={}
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for prefix,folder in sources:
            for path in sorted(folder.rglob('*')):
                if not path.is_file() or '_bin' in path.parts or '__pycache__' in path.parts:continue
                name=prefix+'/'+path.relative_to(folder).as_posix();assert name not in entries
                z.write(path,name);entries[name]=dict(size=path.stat().st_size,sha256=sha(path))
    parts=[]
    with archive.open('rb') as stream:
        while raw:=stream.read(24*1024*1024):
            path=destination/('evidence.zip.part'+str(len(parts)).zfill(3));path.write_bytes(raw)
            parts.append(dict(name=path.name,size=len(raw),sha256=sha(path)))
    index=dict(archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),parts=parts,entries=entries,
               raw_failure_records_retained=True,solver_binaries_archived=False)
    write(destination/'archive-index.json',index)
    return index


def verify(folder):
    m=json.loads((folder/'archive-index.json').read_text(encoding='utf-8'));digest=hashlib.sha256();size=0
    with tempfile.TemporaryFile() as stream:
        for part in m['parts']:
            path=folder/part['name'];assert path.stat().st_size==part['size'] and sha(path)==part['sha256']
            with path.open('rb') as f:
                for raw in iter(lambda:f.read(1024*1024),b''):stream.write(raw);digest.update(raw);size+=len(raw)
        assert size==m['archive_bytes'] and digest.hexdigest()==m['archive_sha256'];stream.seek(0)
        with zipfile.ZipFile(stream) as z:
            assert len(z.namelist())==len(set(z.namelist())) and set(z.namelist())==set(m['entries'])
            for name in z.namelist():
                raw=z.read(name);e=m['entries'][name];assert len(raw)==e['size'] and hashlib.sha256(raw).hexdigest()==e['sha256']
    return dict(parts=len(m['parts']),entries=len(m['entries']),archive_sha256=digest.hexdigest(),all_bytes_verified=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--extract-old',type=Path);p.add_argument('--output',type=Path)
    p.add_argument('--formal',type=Path);p.add_argument('--preparation',type=Path);p.add_argument('--old',type=Path)
    p.add_argument('--verify',type=Path);a=p.parse_args()
    if a.extract_old:print(extract_old(a.extract_old,a.output))
    elif a.verify:print(verify(a.verify))
    else:
        print(seal([('formal',a.formal),('retained-preparation',a.preparation),('old-verified-subset',a.old)],a.output)['archive_bytes'])

if __name__=='__main__':main()
