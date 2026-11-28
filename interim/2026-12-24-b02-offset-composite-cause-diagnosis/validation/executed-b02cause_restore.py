import concurrent.futures, hashlib, json, pathlib, subprocess, time, urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
BASE = '33b185b060c26cf2d0957271391b7bec86e491da'
OUT = ROOT/'work/b02cause-preparation'
OUT.mkdir(exist_ok=True)
rows = []
for line in subprocess.check_output(['git','ls-tree','-r',BASE],cwd=ROOT,text=True).splitlines():
    info, path = line.split('\t',1)
    mode, kind, sha = info.split()
    if kind == 'blob' and any(s in pathlib.Path(path).name for s in ('.zip','.tar','.gz','.npz','.zst','.bin')):
        rows.append((path,sha))

def restore(item):
    path, expected = item
    dest = ROOT/path
    start = time.monotonic()
    if dest.exists():
        data = dest.read_bytes()
        actual = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if actual == expected:
            return dict(path=path,git_blob_sha1=expected,size=len(data),already_present=True)
        raise ValueError('existing file differs: '+path)
    for attempt in range(1,5):
        try:
            url = 'https://raw.githubusercontent.com/tak063495-prog/solidcraft_soil_analysis/'+BASE+'/'+path
            with urllib.request.urlopen(url,timeout=120) as response:
                data = response.read()
            actual = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
            if actual != expected:
                raise ValueError('Git blob mismatch: '+path)
            dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(data)
            return dict(path=path,git_blob_sha1=expected,size=len(data),sha256=hashlib.sha256(data).hexdigest(),already_present=False,attempts=attempt,elapsed_seconds=time.monotonic()-start)
        except Exception as e:
            if attempt==4:
                return dict(path=path,error=repr(e))

started=time.monotonic(); done=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
    for row in pool.map(restore,rows):
        done.append(row)
        if len(done)%100==0:
            print('restored/verified',len(done),'/',len(rows),flush=True)
report=dict(fixed_commit=BASE,archives=len(rows),files=done,errors=[r for r in done if 'error' in r],elapsed_seconds=time.monotonic()-started)
(OUT/'archive-restoration.json').write_text(json.dumps(report,indent=2)+'\n')
print('done',len(done),'errors',len(report['errors']),flush=True)
raise SystemExit(bool(report['errors']))
