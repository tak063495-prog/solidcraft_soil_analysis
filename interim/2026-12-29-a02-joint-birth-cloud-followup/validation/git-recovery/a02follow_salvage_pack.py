"""Recover complete Git objects, leaving the failing suffix separately retained."""
from pathlib import Path
import json,struct,zlib,hashlib,os
p=Path('.git/objects/pack/tmp_pack_zFy4wn');report=json.loads(Path('work/a02follow-failed-pack-checksum.json').read_text());ends=[];types={}
with p.open('rb') as f:
 header=f.read(12);count=struct.unpack('>I',header[8:12])[0]
 for index in range(count):
  start=f.tell()
  try:
   b=f.read(1)
   if not b:raise EOFError('object header')
   value=b[0];kind=(value>>4)&7;size=value&15;shift=4
   while value&128:
    b=f.read(1)
    if not b:raise EOFError('object size')
    value=b[0];size|=(value&127)<<shift;shift+=7
   if kind==6:
    b=f.read(1)
    if not b:raise EOFError('offset delta')
    while b[0]&128:
     b=f.read(1)
     if not b:raise EOFError('offset delta')
   elif kind==7:
    if len(f.read(20))!=20:raise EOFError('ref delta')
   elif kind not in (1,2,3,4):raise ValueError(('unknown type',kind))
   dec=zlib.decompressobj();produced=0
   while not dec.eof:
    chunk=f.read(1024*1024)
    if not chunk:raise EOFError('incomplete zlib object')
    data=dec.decompress(chunk);produced+=len(data)
   f.seek(-len(dec.unused_data),1)
   assert produced==size,(produced,size)
   ends.append(f.tell());types[kind]=types.get(kind,0)+1
   if (index+1)%100==0:print('complete objects',index+1,'offset',f.tell(),flush=True)
  except (EOFError,zlib.error,ValueError) as e:
   report.update(incomplete_object_index=index,incomplete_object_offset=start,error=repr(e));break
 cut=ends[-1] if ends else 12
 f.seek(cut);tail=f.read()
Path('work/a02follow-failed-pack-tail.bin').write_bytes(tail)
report.update(complete_objects=len(ends),complete_object_types=types,retained_suffix_bytes=len(tail),retained_suffix_sha256=hashlib.sha256(tail).hexdigest(),complete_prefix_end=cut)
Path('work/a02follow-pack-salvage-scan.json').write_text(json.dumps(report,indent=2)+'\n')
assert len(ends)>=count-2,report
assert types.get(1,0)==types.get(2,0)==types.get(4,0)==0
p.chmod(0o600)
with p.open('r+b') as f:
 f.seek(8);f.write(struct.pack('>I',len(ends)));f.truncate(cut);f.seek(0);h=hashlib.sha1();left=cut
 while left:
  b=f.read(min(left,8*1024*1024));assert b;h.update(b);left-=len(b)
 f.seek(cut);f.write(h.digest())
p.chmod(0o400)
report.update(recovered_pack_sha1=h.hexdigest(),recovered_bytes=cut+20,header_count_adjusted_only=True,working_files_or_commits_changed=False)
Path('work/a02follow-pack-salvage.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
