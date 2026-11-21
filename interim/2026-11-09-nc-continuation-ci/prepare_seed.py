"""Verify every byte/entry of original archives; retain only restart prerequisites."""
import argparse, concurrent.futures, hashlib, json, pathlib, tarfile, time
HERE = pathlib.Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--archives-root', type=pathlib.Path, default=HERE.parent)
parser.add_argument('--output', type=pathlib.Path, required=True)
args = parser.parse_args()
ROOT = args.archives_root.resolve()
OUT = args.output.resolve()
class Parts:
    def __init__(self, paths): self.paths=iter(paths); self.f=None
    def read(self,n=-1):
        if n<0: raise ValueError('bounded read required')
        result=bytearray()
        while len(result)<n:
            if self.f is None:
                try: self.f=next(self.paths).open('rb')
                except StopIteration: break
            b=self.f.read(n-len(result))
            if b: result.extend(b)
            else: self.f.close(); self.f=None
        return bytes(result)
def verify(name):
    start=time.monotonic(); src=ROOT/name; out=OUT/name; out.mkdir(parents=True,exist_ok=False)
    m=json.loads((src/'archive-sha256.json').read_text()); h=hashlib.sha256()
    for p in m['parts']:
        ph=hashlib.sha256(); length=0
        with (src/p['name']).open('rb') as f:
            while b:=f.read(1024*1024): h.update(b); ph.update(b); length+=len(b)
        assert length==p['bytes'] and ph.hexdigest()==p['sha256']
    assert h.hexdigest()==m['archive_sha256']
    seen=set(); retained={}
    with tarfile.open(fileobj=Parts([src/p['name'] for p in m['parts']]),mode='r|xz') as z:
        for ent in z:
            p=pathlib.PurePosixPath(ent.name)
            assert ent.isfile() and not p.is_absolute() and '..' not in p.parts and ent.name not in seen
            keep=(name.endswith('further-branches') and (ent.name.startswith('observer-source/') or ent.name in ['extended-native-further/restart.json','extended-independent-further/restart.json','extended-native-further/controls.json','extended-independent-further/controls.json'])) or (name.endswith('beyond-fold') and ent.name.startswith('origin/')) or (name.endswith('pseudo-arclength') and ent.name=='independent-arc-history400/samples.json')
            target=out/ent.name
            if keep: target.parent.mkdir(parents=True,exist_ok=True); dest=target.open('xb')
            else: dest=None
            eh=hashlib.sha256()
            with z.extractfile(ent) as f:
                while b:=f.read(1024*1024):
                    eh.update(b)
                    if dest: dest.write(b)
            if dest: dest.close(); retained[ent.name]=eh.hexdigest()
            assert eh.hexdigest()==m['entries'][ent.name],ent.name
            seen.add(ent.name)
    assert seen==set(m['entries'])
    r=dict(passed=False,solver_executed=False,archive_sha256=h.hexdigest(),verified_entries=len(seen),retained=retained,elapsed_seconds=time.monotonic()-start)
    (out/'verification.json').write_text(json.dumps(r,indent=2)+'\n')
    print(name,len(seen),r['elapsed_seconds'],flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    list(pool.map(verify,['2026-11-09-nc-further-branches','2026-11-09-nc-beyond-fold','2026-11-09-nc-pseudo-arclength']))
