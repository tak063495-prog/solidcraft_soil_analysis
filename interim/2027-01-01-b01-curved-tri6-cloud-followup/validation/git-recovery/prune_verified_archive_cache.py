"""Remove only refetchable single-archive promisor packs after byte verification."""
from pathlib import Path
import hashlib,json,subprocess,shutil
out=Path('work/b01follow-git-recovery');available=json.loads((out/'fetch-result-v3.json').read_text())
assert available['all_required_objects_available'] and not available['remaining']
materialized=json.loads((out/'materialization.json').read_text())['files'];verified={}
for row in materialized:
 p=Path(row['path']);h=hashlib.sha1(('blob '+str(p.stat().st_size)+'\0').encode());s=hashlib.sha256()
 with p.open('rb') as f:
  for raw in iter(lambda:f.read(1024*1024),b''):h.update(raw);s.update(raw)
 assert h.hexdigest()==row['git_blob_sha1'] and s.hexdigest()==row['sha256'],row['path']
 verified.setdefault(h.hexdigest(),[]).append(row)
removed=[];before=shutil.disk_usage('.').free
for idx in sorted(Path('.git/objects/pack').glob('pack-*.idx')):
 pack=idx.with_suffix('.pack');promisor=idx.with_suffix('.promisor')
 if not pack.exists() or not promisor.exists() or pack.stat().st_size<1024*1024:continue
 with idx.open('rb') as f:lines=subprocess.check_output(['git','show-index'],stdin=f,text=True).splitlines()
 if len(lines)!=1:continue
 oid=lines[0].split()[1]
 if oid not in verified:continue
 with pack.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 removed.append(dict(pack=str(pack),bytes=pack.stat().st_size,pack_sha256=digest,git_blob_sha1=oid,verified_work_files=verified[oid]))
 pack.unlink();idx.unlink();promisor.unlink();idx.with_suffix('.rev').unlink(missing_ok=True)
(out/'cache-pruning.json').write_text(json.dumps(dict(only_refetchable_single_blob_promisor_cache=True,working_files_or_commits_deleted=False,removed=removed,
 freed_bytes=shutil.disk_usage('.').free-before,free_bytes=shutil.disk_usage('.').free),indent=2)+'\n')
print('pruned',len(removed),'free',shutil.disk_usage('.').free,flush=True)
