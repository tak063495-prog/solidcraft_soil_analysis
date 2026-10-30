"""Verify the actual published635 CI artifact against every source and lock."""
from pathlib import Path,PurePosixPath
import hashlib,json,re,zipfile
p=Path(__file__).resolve().parent
root=p.parents[1]
out=p/'Root635-remote-CI-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
commit='458c260da9ad286332cac60e50dfc510245d0cde'
meta=read(out/'remote-metadata.json')
assert meta['run']['id']==40983704233 and meta['job']['id']==123829597537
assert meta['run']['status']=='completed' and meta['run']['conclusion']=='success'
assert meta['job']['status']=='completed' and meta['job']['conclusion']=='success'
assert meta['run']['head_sha']==meta['artifact']['head_sha']==commit
assert all(s['conclusion']=='success' for s in meta['job']['steps'])
archive=out/'artifact.zip';digest=sha(archive)
assert archive.stat().st_size==meta['artifact']['bytes']
assert 'sha256:'+digest==meta['artifact']['digest']
publication=root/'runs/interim-upload-20261030-v49'
receipt=read(publication/'publication-receipt.json')
plan=read(publication/'upload-plan.json')
assert receipt['passed'] and receipt['commit_sha']==commit
names={'python-packages.log','rust-tests.log','python-tests.log','rust-version.log','source-snapshot.json'}
with zipfile.ZipFile(archive) as z:
    assert len(z.namelist())==len(set(z.namelist()))==5 and set(z.namelist())==names
    assert z.testzip() is None and sum(i.file_size for i in z.infolist())<2*(1<<20)
    files={name:z.read(name) for name in names}
    for name,raw in files.items():
        q=PurePosixPath(name);assert not q.is_absolute() and '..' not in q.parts
        with (out/name).open('xb') as f:f.write(raw)
source=json.loads(files['source-snapshot.json'].decode('utf-8-sig'))
frozen=read(p/'Root-source635-promotion-v1/quality-v1/source-snapshot.json')
assert len(source)==635 and source.keys()==frozen.keys()
exact=[];crlf=[];local={}
for name,h in source.items():
    q=publication/'publication'/name
    if not q.exists():q=root/name
    raw=q.read_bytes()
    expected=plan['files'][name]['sha256']
    assert hashlib.sha256(raw).hexdigest()==expected,'published reference changed: '+name
    local[name]=expected
    if expected==h:exact.append(name)
    else:
        assert b'\r' not in raw.replace(b'\r\n',b''),'standalone CR in published source'
        windows=raw.replace(b'\r\n',b'\n').replace(b'\n',b'\r\n')
        assert hashlib.sha256(windows).hexdigest()==h,'unexplained CI source difference: '+name
        crlf.append(name)
dump(out/'source-checkout-line-endings.json',dict(passed=True,files=635,
    byte_exact=exact,exact_CRLF_transform=crlf,actual_CI_manifest_retained=True,
    published_Root_source_hashes=local,no_arbitrary_whitespace_normalization=True,
    published_docs_updated_after_local219_test_source_archive=True))
rust=files['rust-tests.log'].decode('utf-8-sig').replace('\r','')
summaries=[tuple(map(int,row)) for row in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored',rust)]
assert len(summaries)==48 and all(failed==ignored==0 for _,failed,ignored in summaries)
assert not re.search(r'test result: FAILED|error: test failed',rust)
assert sum(s[0] for s in summaries)==300
python=files['python-tests.log'].decode('utf-8-sig').replace('\r','')
counts=re.findall(r'^Ran (\d+) tests? in ',python,re.MULTILINE)
assert counts==['219'] and re.search(r'^OK\s*$',python,re.MULTILINE)
assert not re.search(r'^FAILED\b|\.\.\. skipped\b|^OK \(skipped=',python,re.MULTILINE)
packages=dict(line.strip().split('==',1) for line in files['python-packages.log'].decode('utf-8-sig').splitlines() if '==' in line)
required={}
for name in ['requirements.lock.txt','requirements-reference.lock.txt']:
    q=publication/'publication'/name
    if not q.exists():q=root/name
    assert sha(q)==plan['files'][name]['sha256']
    for line in q.read_text('utf-8-sig').splitlines():
        if '==' in line:
            key,value=line.strip().split('==',1);required[key]=value
assert packages==required and len(packages)==9 and packages['ijson']=='3.6.0'
version=files['rust-version.log'].decode('utf-8-sig')
assert 'rustc 1.99.1' in version
dump(out/'terminal.json',dict(passed=True,commit_sha=commit,run_id=meta['run']['id'],job_id=meta['job']['id'],
    artifact_id=meta['artifact']['id'],artifact_bytes=archive.stat().st_size,artifact_sha256=digest,
    all5_archive_members_CRC_hash_and_paths_verified=True,all635_published_sources_verified=True,
    byte_exact_files=len(exact),exact_checkout_CRLF_files=len(crlf),byte_exact_checkout_not_claimed=True,
    Rust_tests_passed=300,Rust_tests_failed=0,Rust_tests_ignored=0,Rust_summary_count=48,
    Python_tests_passed=219,Python_tests_skipped=0,all9_pinned_dependencies=packages,
    all_remote_steps_success=True,actual_remote_published635_not_prior621=True,
    current_Native_async_UI_or_new_PC_other_GPU_distribution_or_general_convergence_not_proven=True,
    p6_complete=False))
print('PASS actual remote published Root635: Rust300 Python219 dependencies9 and all635 source identities',flush=True)
