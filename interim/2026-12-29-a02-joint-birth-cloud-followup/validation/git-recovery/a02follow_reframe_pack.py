from pathlib import Path
import json,hashlib,struct
r=json.loads(Path('work/a02follow-pack-salvage-scan.json').read_text());p=Path(r['path']);assert p.stat().st_size==r['size'];assert set(map(int,r['complete_object_types']))<={3,6};assert r['complete_objects']>0
cut=r['complete_prefix_end'];old_header=bytes.fromhex(r['header_hex']);new_header=old_header[:8]+struct.pack('>I',r['complete_objects']);new_hash=hashlib.sha1(new_header);old_hash=hashlib.sha256(old_header)
with p.open('rb') as f:
 assert f.read(12)==old_header;left=cut-12
 while left:
  b=f.read(min(left,8*1024*1024));assert b;new_hash.update(b);old_hash.update(b);left-=len(b)
 suffix=f.read();assert hashlib.sha256(suffix).hexdigest()==r['retained_suffix_sha256'];old_hash.update(suffix)
p.chmod(0o600)
with p.open('r+b') as f:f.write(new_header);f.truncate(cut);f.seek(cut);f.write(new_hash.digest())
p.chmod(0o400);target=p.with_name('pack-'+new_hash.hexdigest()+'.pack');p.rename(target)
r.update(original_failed_pack_sha256=old_hash.hexdigest(),recovered_pack_sha1=new_hash.hexdigest(),recovered_pack=str(target),recovered_bytes=cut+20,working_files_or_commits_changed=False,retained_original_header_and_incomplete_suffix=True)
Path('work/a02follow-pack-salvage.json').write_text(json.dumps(r,indent=2)+'\n');print(str(target),flush=True)
