from pathlib import Path
import hashlib,struct,json
p=Path('.git/objects/pack/tmp_pack_zFy4wn');size=p.stat().st_size;h=hashlib.sha1();left=size-20
with p.open('rb') as f:
 header=f.read(12);f.seek(0)
 while left:
  b=f.read(min(left,8*1024*1024));assert b;h.update(b);left-=len(b)
 tail=f.read()
r=dict(path=str(p),size=size,header_hex=header.hex(),object_count=struct.unpack('>I',header[8:12])[0],expected_sha1=h.hexdigest(),trailer_hex=tail.hex(),valid_complete_pack=h.digest()==tail)
Path('work/a02follow-failed-pack-checksum.json').write_text(json.dumps(r,indent=2)+'\n');print(r)
