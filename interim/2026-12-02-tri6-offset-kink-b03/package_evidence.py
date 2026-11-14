"""Lossless split evidence with every part/member hash; verification is not a solve."""
import argparse,hashlib,json,shutil,zipfile
from pathlib import Path
from b03_cases import HERE,ROOT,read,write


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def create(work,out):
    archive=work/'b03-evidence.zip';members=[]
    groups={'fresh':'b03-fresh','probes':'b03-probes-current','supplemental':'b03-supplemental',
            'mutations':'b03-mutations','retained-initial-mutations':'b03-mutations-initial-error','retained-prebuild-probes':'b03-probes'}
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for label,folder in groups.items():
            for p in sorted((work/folder).rglob('*')):
                if not p.is_file() or '__pycache__' in p.parts:continue
                arc=label+'/'+str(p.relative_to(work/folder)).replace('\\','/')
                z.write(p,arc);members.append(dict(path=arc,bytes=p.stat().st_size,sha256=sha(p)))
        for p in sorted(work.glob('pr-*.txt'))+sorted(work.glob('pr-*.html')):
            arc='reference-pr-pages/'+p.name;z.write(p,arc);members.append(dict(path=arc,bytes=p.stat().st_size,sha256=sha(p)))
    parts=[]
    with archive.open('rb') as f:
        index=0
        while data:=f.read(16*1024*1024):
            p=out/f'evidence.zip.part{index:03d}'
            with p.open('xb') as target:target.write(data)
            parts.append(dict(path=p.name,bytes=len(data),sha256=sha(p)));index+=1
    manifest=dict(schema='b03-split-evidence/1',parts=parts,archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),members=members,
         no_recomputed_measurement_is_a_fresh_solve=True)
    write(out/'evidence-manifest.json',manifest)
    binding=read(work/'b03-fresh/source-binding-before.json');snapshot=out/'source-snapshot.zip';source_members=[]
    paths=set(binding)|{str(p.relative_to(ROOT)) for p in HERE.rglob('*') if p.is_file() and p.suffix in ('.py','.rs','.md')}
    paths|={str(p.relative_to(ROOT)) for p in ROOT.glob('requirements*.txt')}
    paths|={str(p.relative_to(ROOT)) for p in (ROOT/'.github/workflows').glob('*')}
    paths.add('sources.lock.json')
    with zipfile.ZipFile(snapshot,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for relative in sorted(paths):
            p=ROOT/relative;digest=sha(p)
            if relative in binding:assert digest==binding[relative]
            z.write(p,relative);source_members.append(dict(path=relative,sha256=digest,bytes=p.stat().st_size,at_start_of_fresh_solve=relative in binding))
    write(out/'source-snapshot-manifest.json',dict(sha256=sha(snapshot),bytes=snapshot.stat().st_size,members=source_members,
            fresh_solve_sources_bound_before_and_after=True,postprocessing_and_prose_added_after_start_separately_labeled=True))
    return manifest


def verify(manifest_path,scratch):
    d=read(manifest_path);scratch.mkdir(parents=True,exist_ok=False);archive=scratch/'combined.zip'
    with archive.open('xb') as target:
        for part in d['parts']:
            p=manifest_path.parent/part['path'];assert p.stat().st_size==part['bytes'] and sha(p)==part['sha256']
            with p.open('rb') as f:shutil.copyfileobj(f,target)
    assert archive.stat().st_size==d['archive_bytes'] and sha(archive)==d['archive_sha256']
    members={v['path']:v for v in d['members']}
    with zipfile.ZipFile(archive) as z:
        assert set(z.namelist())==members.keys() and len(z.namelist())==len(members)
        for name,v in members.items():
            h=hashlib.sha256();size=0
            with z.open(name) as f:
                for b in iter(lambda:f.read(1024*1024),b''):h.update(b);size+=len(b)
            assert size==v['bytes'] and h.hexdigest()==v['sha256'],name
    result=dict(passed=True,archive_sha256=d['archive_sha256'],parts_verified=len(d['parts']),members_verified=len(members),
                verification_is_not_a_new_solve=True)
    write(scratch/'verification.json',result);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path);p.add_argument('--output',type=Path,default=HERE)
    p.add_argument('--verify',type=Path);p.add_argument('--scratch',type=Path);a=p.parse_args()
    if a.verify:print(json.dumps(verify(a.verify,a.scratch)))
    else:create(a.work,a.output)
