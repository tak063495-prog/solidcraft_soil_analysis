import subprocess,pathlib,json,time,hashlib
root=pathlib.Path('/workspace/solidcraft_soil_analysis');base='7b254b1cdbe9c288989631944cd698aa861f83ac'
rows=[]
for line in subprocess.check_output(['git','ls-tree','-r',base],cwd=root,text=True).splitlines():
 meta,path=line.split('\t');mode,kind,oid=meta.split()
 if kind=='blob' and not (root/path).exists():rows.append((oid,path))
print('Missing original files',len(rows),flush=True)
for i in range(0,len(rows),100):
 batch=rows[i:i+100];oids=sorted(set(x[0] for x in batch))
 p=subprocess.run(['git','-c','fetch.negotiationAlgorithm=noop','fetch','--no-tags','--no-write-fetch-head','--no-auto-gc','origin','--stdin'],cwd=root,input='\n'.join(oids)+'\n',text=True,stdout=subprocess.DEVNULL)
 if p.returncode:raise RuntimeError(p.returncode)
 for oid,name in batch:
  path=root/name;path.parent.mkdir(parents=True,exist_ok=True)
  with path.open('xb') as f:subprocess.run(['git','cat-file','blob',oid],cwd=root,stdout=f,check=True)
  raw=path.read_bytes();actual=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest();assert actual==oid,name
 print('Restored',min(i+100,len(rows)),flush=True)
(root/'work/a02-materialization.json').write_text(json.dumps({'base':base,'restored':rows,'all_restored_blobs_verified':True},indent=2)+'\n')
