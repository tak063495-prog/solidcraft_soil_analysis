import hashlib,json,pathlib,subprocess,time,concurrent.futures
root=pathlib.Path.cwd();base='3747630ab615b4be2117ca1aaa13e798cd436374';items=[]
for line in subprocess.check_output(['git','ls-tree','-r',base],text=True).splitlines():
 info,name=line.split('\t',1);mode,kind,expected=info.split()
 if kind=='blob':items.append((name,mode,expected))
def check(item):
 name,mode,expected=item;p=root/name
 if not p.is_file():return dict(path=name,missing=True)
 size=p.stat().st_size;h=hashlib.sha1(b'blob '+str(size).encode()+b'\0');s=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b);s.update(b)
 return dict(path=name,size=size,expected_git_blob=expected,actual_git_blob=h.hexdigest(),sha256=s.hexdigest(),same=h.hexdigest()==expected,mode_matches=bool(p.stat().st_mode&0o111)==(mode=='100755'))
t=time.monotonic()
with concurrent.futures.ThreadPoolExecutor(4) as pool:rows=list(pool.map(check,items))
report=dict(base=base,files=len(rows),bytes=sum(r.get('size',0) for r in rows),missing=[r['path'] for r in rows if r.get('missing')],differences=[r['path'] for r in rows if not r.get('missing') and not r.get('same')],mode_differences=[r['path'] for r in rows if not r.get('missing') and not r.get('mode_matches')],elapsed_seconds=time.monotonic()-t,rows=rows)
path=root/'work/b03cause-after-preservation.json';path.write_text(json.dumps(report,indent=2)+'\n');print({k:v for k,v in report.items() if k!='rows'},flush=True)
