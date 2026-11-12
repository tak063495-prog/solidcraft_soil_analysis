import concurrent.futures,glob,hashlib,json,os,struct,subprocess,time,urllib.parse,zlib
from pathlib import Path
ROOT=Path.cwd();START='3b1133d3be9905dc24f62c27b5b1eb456a83f2e2'
out=ROOT/'work/artifact-preparation';out.mkdir(exist_ok=True)
tree={}
for line in subprocess.check_output(['git','ls-tree','-r','HEAD'],text=True).splitlines():
 meta,name=line.split('\t');mode,kind,oid=meta.split()
 if kind=='blob':tree[name]=oid
wanted={}
for label in ('2026-11-10-nc-mean-stress-release','2026-11-10-nc-spatial-release','2026-11-15-wet-joint-convergence','2026-11-18-wet-joint-step-convergence'):
 folder=ROOT/'interim'/label;m=json.loads((folder/'archive-sha256.json').read_text())
 for part in m['parts']:
  path=(folder/part['name']).relative_to(ROOT).as_posix()
  wanted[path]=dict(git_blob=tree[path],sha256=part['sha256'],size=part.get('size',part.get('bytes')))
byid={v['git_blob']:k for k,v in wanted.items()};records=[]
def inspect(path,expected):
 h=hashlib.sha256();g=hashlib.sha1(f'blob {path.stat().st_size}\0'.encode())
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b);g.update(b)
 assert path.stat().st_size==expected['size'] and h.hexdigest()==expected['sha256'] and g.hexdigest()==expected['git_blob'],str(path)
for path,v in wanted.items():
 p=ROOT/path
 if p.is_file():inspect(p,v);records.append(dict(path=path,source='already placed exact original Git blob',**v))
# Salvage only complete non-delta blobs from the interrupted network pack;
# neither reference values nor archived snapshots are synthesized.
for pack in glob.glob(str(ROOT/'.git/objects/pack/tmp_pack_*')):
 completed=0;found=0
 try:
  with open(pack,'rb') as f:
   magic,version,count=struct.unpack('>4sII',f.read(12));assert magic==b'PACK' and version in (2,3)
   for number in range(count):
    b=f.read(1)
    if not b:raise EOFError('unfinished object header')
    b=b[0];kind=(b>>4)&7;size=b&15;shift=4
    while b&128:
     b=f.read(1)[0];size|=(b&127)<<shift;shift+=7
    if kind==6:
     b=f.read(1)[0]
     while b&128:b=f.read(1)[0]
    elif kind==7:f.read(20)
    h=hashlib.sha1(f'blob {size}\0'.encode());payload=out/'salvage.tmp';store=kind==3 and size>=1024*1024
    z=zlib.decompressobj();stream=payload.open('wb') if store else None;length=0
    try:
     while not z.eof:
      block=f.read(65536)
      if not block:raise EOFError('unfinished compressed object')
      pending=block
      while pending:
       data=z.decompress(pending,1024*1024);pending=z.unconsumed_tail
       if store:stream.write(data);h.update(data);length+=len(data)
       if z.eof:break
      if z.eof:f.seek(-len(z.unused_data),os.SEEK_CUR)
    finally:
     if stream:stream.close()
    completed+=1
    if store:
     oid=h.hexdigest()
     if oid in byid:
      path=byid[oid];target=ROOT/path
      if not target.exists():
       inspect(payload,wanted[path]);target.parent.mkdir(parents=True,exist_ok=True);payload.rename(target);found+=1
       records.append(dict(path=path,source='complete blob salvaged from interrupted fetch; SHA256 and original Git blob verified',**wanted[path]))
     if payload.exists():payload.unlink()
 except Exception as e:print('salvage ended',type(e).__name__,str(e),flush=True)
 print('salvage',completed,'objects',found,'needed parts',flush=True)
 if (out/'salvage.tmp').exists():(out/'salvage.tmp').unlink()
(out/'salvage.json').write_text(json.dumps(records,indent=2)+'\n')
missing=[p for p in wanted if not (ROOT/p).is_file()]
print('need download',len(missing),'parts',sum(wanted[p]['size'] for p in missing),'bytes',flush=True)
def get(path):
 v=wanted[path];target=ROOT/path;temp=out/(v['git_blob']+'.download');log=out/(v['git_blob']+'.log')
 url='https://raw.githubusercontent.com/tak063495-prog/solidcraft_soil_analysis/'+START+'/'+urllib.parse.quote(path)
 start=time.monotonic()
 with log.open('x') as stream:
  p=subprocess.run(['curl','-fSL','--retry','2','--retry-delay','2','--connect-timeout','30','--max-time','240',url,'-o',str(temp)],stdout=stream,stderr=subprocess.STDOUT)
 assert p.returncode==0,(path,p.returncode)
 inspect(temp,v);assert not target.exists();target.parent.mkdir(parents=True,exist_ok=True);temp.rename(target)
 return dict(path=path,source='raw GitHub at fixed commit, SHA256 and original Git blob verified',elapsed_seconds=time.monotonic()-start,**v)
with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
 for n,future in enumerate(concurrent.futures.as_completed([pool.submit(get,p) for p in missing])):
  try:records.append(future.result())
  except Exception as e:records.append(dict(error=repr(e)))
  if n%16==0:print('download completed',n+1,'/',len(missing),flush=True)
  (out/'results.json').write_text(json.dumps(records,indent=2)+'\n')
assert all((ROOT/p).is_file() for p in wanted) and not any('error' in r for r in records)
print('all original artifact parts placed and verified',len(wanted),flush=True)
