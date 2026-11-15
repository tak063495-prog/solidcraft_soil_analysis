from pathlib import Path
import json,struct,zlib,hashlib,subprocess
out=Path('work/b01follow-git-recovery');p=Path('.git/objects/pack/tmp_pack_t0zvNv');ends=[];types={};before_size=p.stat().st_size
with p.open('rb') as f:
 header=f.read(12);count=struct.unpack('>I',header[8:12])[0]
 for index in range(count):
  start=f.tell()
  try:
   b=f.read(1)
   if not b:raise EOFError('header')
   v=b[0];kind=(v>>4)&7;size=v&15;shift=4
   while v&128:
    b=f.read(1)
    if not b:raise EOFError('size')
    v=b[0];size|=(v&127)<<shift;shift+=7
   if kind==6:
    b=f.read(1)
    if not b:raise EOFError('offset')
    while b[0]&128:
     b=f.read(1)
     if not b:raise EOFError('offset')
   elif kind==7:
    if len(f.read(20))!=20:raise EOFError('ref')
   elif kind not in (1,2,3,4):raise ValueError(('type',kind))
   dec=zlib.decompressobj();produced=0
   while not dec.eof:
    b=f.read(1024*1024)
    if not b:raise EOFError('incomplete zlib object')
    produced+=len(dec.decompress(b))
   f.seek(-len(dec.unused_data),1)
   if produced!=size:raise ValueError(('size',produced,size))
   ends.append(f.tell());types[kind]=types.get(kind,0)+1
   if (index+1)%100==0:print('complete',index+1,'offset',f.tell(),flush=True)
  except (EOFError,zlib.error,ValueError) as e:
   error=repr(e);break
 cut=ends[-1] if ends else 12;f.seek(cut);tail=f.read()
(out/'original-header.bin').write_bytes(header);(out/'incomplete-suffix.bin').write_bytes(tail)
r=dict(original_path=str(p),original_size=before_size,original_header_hex=header.hex(),advertised_objects=count,
 complete_objects=len(ends),complete_types=types,cut=cut,error=error,suffix_bytes=len(tail),suffix_sha256=hashlib.sha256(tail).hexdigest())
assert ends and set(types)<={3,6},r
newheader=header[:8]+struct.pack('>I',len(ends));h=hashlib.sha1(newheader);original_sha=hashlib.sha256(header)
with p.open('rb') as f:
 f.seek(12);left=cut-12
 while left:
  b=f.read(min(left,8*1024*1024));assert b;left-=len(b);h.update(b);original_sha.update(b)
 original_sha.update(tail)
p.chmod(0o600)
with p.open('r+b') as f:f.write(newheader);f.truncate(cut);f.seek(cut);f.write(h.digest())
newpath=p.with_name('pack-'+h.hexdigest()+'.pack');p.rename(newpath)
r.update(original_sha256=original_sha.hexdigest(),recovered_pack=str(newpath),recovered_pack_sha1=h.hexdigest(),only_header_count_and_incomplete_suffix_changed=True,working_files_or_commits_changed=False)
(out/'recovery.json').write_text(json.dumps(r,indent=2)+'\n')
with (out/'index-pack.log').open('x') as f:
 subprocess.run(['git','index-pack','--strict','--threads=2',str(newpath)],stdout=f,stderr=subprocess.STDOUT,check=True)
newpath.with_suffix('.promisor').touch()
print('recovered and strict indexed',len(ends),'objects',flush=True)
