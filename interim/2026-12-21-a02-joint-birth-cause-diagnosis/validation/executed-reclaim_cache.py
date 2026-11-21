import pathlib,subprocess,json,time,shutil
root=pathlib.Path('/workspace/solidcraft_soil_analysis')
while not (root/'work/a02-materialization.json').exists():time.sleep(1)
known=json.loads((root/'work/a02-materialization.json').read_text())['restored'];oids={o for o,p in known}
removed=[]
for pack in sorted((root/'.git/objects/pack').glob('*.pack'),key=lambda p:p.stat().st_size,reverse=True):
 if shutil.disk_usage(root).free >= 20*1024**3:break
 idx=pack.with_suffix('.idx')
 if not pack.with_suffix('.promisor').exists() or not idx.exists():continue
 with idx.open('rb') as f:listing=subprocess.check_output(['git','show-index'],cwd=root,stdin=f,text=True)
 objects={line.split()[1] for line in listing.splitlines()}
 if not objects or not objects<=oids:continue
 paths=[p for o,p in known if o in objects]
 assert all((root/p).is_file() for p in paths)
 size=sum(p.stat().st_size for p in [pack,idx,pack.with_suffix('.promisor'),pack.with_suffix('.rev')] if p.exists())
 for p in [pack,idx,pack.with_suffix('.promisor'),pack.with_suffix('.rev')]:p.unlink(missing_ok=True)
 removed.append(dict(pack=pack.name,object_count=len(objects),bytes=size,all_objects_are_verified_restored_original_blobs=True))
result=dict(removed_promisor_download_cache_only=removed,tracked_evidence_deleted=False,commits_and_trees_deleted=False,free_bytes=shutil.disk_usage(root).free)
(root/'work/a02-validation/cache-reclamation.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
