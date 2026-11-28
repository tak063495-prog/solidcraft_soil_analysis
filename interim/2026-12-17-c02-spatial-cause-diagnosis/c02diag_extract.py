"""Select saved contact-field chunks after verifying the complete old ZIP seal.

Read split ZIP parts directly; never rewrite the old archive or extract its
large JSONL equilibrium histories merely to measure existing saved fields.
"""
import argparse
import bisect
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OLD=ROOT/'interim/2026-11-19-wet-joint-spatial-convergence'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


class SplitZip(io.RawIOBase):
    def __init__(self,parts):
        self.files=[p.open('rb') for p in parts];self.ends=[];n=0
        for p in parts:n+=p.stat().st_size;self.ends.append(n)
        self.position=0;self.length=n
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        p=offset if whence==0 else self.position+offset if whence==1 else self.length+offset
        if p<0:raise ValueError('negative ZIP seek')
        self.position=p;return p
    def read(self,size=-1):
        remaining=max(0,self.length-self.position) if size<0 else min(size,max(0,self.length-self.position))
        chunks=[]
        while remaining:
            i=bisect.bisect_right(self.ends,self.position);start=0 if i==0 else self.ends[i-1]
            count=min(remaining,self.ends[i]-self.position)
            self.files[i].seek(self.position-start);chunk=self.files[i].read(count)
            if len(chunk)!=count:raise OSError('truncated ZIP part')
            chunks.append(chunk);remaining-=count;self.position+=count
        return b''.join(chunks)
    def close(self):
        for f in self.files:f.close()
        super().close()


def extract(out):
    out.mkdir(parents=True,exist_ok=False)
    index=json.loads((OLD/'archive-sha256.json').read_text(encoding='utf-8'))
    archive_hash=hashlib.sha256();parts=[]
    for item in index['parts']:
        p=OLD/item['file'];h=hashlib.sha256()
        if p.stat().st_size!=item['bytes']:raise ValueError('ZIP part length mismatch')
        with p.open('rb') as source:
            for chunk in iter(lambda:source.read(1024*1024),b''):h.update(chunk);archive_hash.update(chunk)
        if h.hexdigest()!=item['sha256']:raise ValueError('ZIP part SHA mismatch')
        parts.append(p)
    if archive_hash.hexdigest()!=index['zip_sha256']:raise ValueError('complete ZIP SHA mismatch')
    selected={}
    def wanted(name):
        if name in {f'paths/{n}.json' for n in ('mesh-h0','mesh-h1','mesh-h2','ground','stages','identity')}:return True
        match=re.fullmatch(r'paths/h([012])-g(32|64|128|256)/(native|independent)/(manifest.json|[0-9]+.npz)',name)
        return bool(match and not (match[1]=='2' and (match[2] not in ('32','64') or match[3]=='independent')))
    with SplitZip(parts) as raw,zipfile.ZipFile(raw) as archive:
        if set(archive.namelist())!=set(index['entries']):raise ValueError('ZIP entry inventory mismatch')
        for name in sorted(index['entries']):
            if not wanted(name):continue
            if '..' in Path(name).parts:raise ValueError('unsafe ZIP name')
            target=out/Path(name).relative_to('paths');target.parent.mkdir(parents=True,exist_ok=True)
            with archive.open(name) as source,target.open('xb') as stream:
                while chunk:=source.read(1024*1024):stream.write(chunk)
            expected=index['entries'][name]
            if target.stat().st_size!=expected['bytes'] or digest(target)!=expected['sha256']:raise ValueError('selected entry SHA mismatch')
            selected[name]=expected
    result=dict(old_archive_index_sha256=digest(OLD/'archive-sha256.json'),old_zip_sha256=index['zip_sha256'],all_parts_verified=len(parts),complete_zip_verified=True,selected_entries=selected,excluded_raw_histories_not_used_as_new_solves=True)
    with (out/'selection-provenance.json').open('x',encoding='utf-8',newline='\n') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(selected_entries=len(selected),selected_bytes=sum(x['bytes'] for x in selected.values()),old_zip_sha256=index['zip_sha256'])),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();extract(a.output.resolve())
