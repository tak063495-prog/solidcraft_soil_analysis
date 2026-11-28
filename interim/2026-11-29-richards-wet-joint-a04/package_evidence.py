"""Lossless split ZIP + every-entry SHA; verification never counts as new solve."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile
HERE=Path(__file__).resolve().parent


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path);p.add_argument('--destination',type=Path,default=HERE);p.add_argument('--labels',nargs='+');p.add_argument('--verify',type=Path);p.add_argument('--scratch',type=Path,required=True);a=p.parse_args();scratch=a.scratch.resolve();scratch.mkdir(parents=True,exist_ok=False)
    if a.verify:
        manifest=json.loads(a.verify.read_text());here=a.verify.parent;archive=scratch/'evidence.zip'
        with archive.open('xb') as out:
            for part in manifest['parts']:
                path=here/part['name'];assert path.stat().st_size==part['bytes'] and digest(path)==part['sha256']
                with path.open('rb') as src:shutil.copyfileobj(src,out)
        assert archive.stat().st_size==manifest['zip_bytes'] and digest(archive)==manifest['zip_sha256']
    else:
        work=a.work.resolve();selected=[]
        labels=a.labels or ['formal-soil-03','formal-soil-02','formal-dry-01','formal-hydrostatic-01','formal-refusals-01','development-01','development-02','development-03','development-refusals-01','development-dry-01','formal-soil-01','numerical-derivative-diagnostic']
        for label in labels:
            for path in (work/label).rglob('*'):
                if path.is_file() and 'bin' not in path.relative_to(work).parts and '__pycache__' not in path.parts:selected.append((path,'runs/'+path.relative_to(work).as_posix()))
        for path in work.iterdir():
            if path.is_file() and path.suffix in ('.log','.json','.py'):selected.append((path,'preparation/'+path.name))
        path=work/'artifacts/results.json'
        if path.is_file():selected.append((path,'acquisition/original-blob-restoration.json'))
        for path in HERE.iterdir():
            if path.is_file() and (path.suffix=='.py' or path.name=='SPEC.ja.md'):selected.append((path,'supplemental-sources/'+path.name))
        entries={};archive=scratch/'evidence.zip'
        with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for path,name in sorted(selected,key=lambda x:x[1]):
                assert name not in entries
                entries[name]=dict(bytes=path.stat().st_size,sha256=digest(path));z.write(path,name)
        parts=[]
        with archive.open('rb') as f:
            for index in range(100000):
                b=f.read(24*1024*1024)
                if not b:break
                name=f'evidence.zip.part{index:03}';a.destination.mkdir(parents=True,exist_ok=True);dest=a.destination/name
                with dest.open('xb') as out:out.write(b)
                parts.append(dict(name=name,bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
        manifest=dict(zip_bytes=archive.stat().st_size,zip_sha256=digest(archive),parts=parts,entries=entries,binaries_excluded_but_SHA_recorded=True,saved_reaudit_not_fresh_solve=True)
        with (a.destination/'evidence-manifest.json').open('x') as f:json.dump(manifest,f,ensure_ascii=False,indent=2);f.write('\n')
    with zipfile.ZipFile(archive) as z:
        assert set(z.namelist())==set(manifest['entries']) and len(z.namelist())==len(set(z.namelist()))
        for info in z.infolist():
            assert not info.is_dir() and not Path(info.filename).is_absolute() and '..' not in Path(info.filename).parts
            h=hashlib.sha256();count=0
            with z.open(info) as stream:
                for b in iter(lambda:stream.read(1024*1024),b''):h.update(b);count+=len(b)
            expected=manifest['entries'][info.filename];assert count==expected['bytes'] and h.hexdigest()==expected['sha256'],info.filename
    result=dict(complete=True,parts=len(manifest['parts']),entries=len(manifest['entries']),bytes=manifest['zip_bytes'],every_entry_sha_verified=True,new_solve_count=0)
    with (scratch/'verification.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
