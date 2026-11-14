"""Verify all new evidence, then remeasure with sealed sources in isolated UTF-8."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

HERE=Path(__file__).resolve().parent


def sha(raw):return hashlib.sha256(raw).hexdigest()


def extract_checked(archive,destination,expected):
    seen=set()
    with zipfile.ZipFile(archive) as z:
        for item in z.infolist():
            name=item.filename;p=Path(name)
            if p.is_absolute() or '..' in p.parts or item.is_dir() or name in seen or name not in expected:
                raise ValueError('unsafe/unlisted/duplicate evidence member')
            raw=z.read(item);spec=expected[name]
            if sha(raw)!=spec['sha256'] or len(raw)!=spec['size']:raise ValueError('changed evidence member')
            target=destination/p;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw);seen.add(name)
    if seen!=set(expected):raise ValueError('missing evidence member')


def run(output):
    output.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((HERE/'archive-sha256.json').read_text(encoding='utf-8'))
    archive=output/'joined.zip';digest=hashlib.sha256();size=0
    with archive.open('xb') as stream:
        for item in manifest['parts']:
            raw=(HERE/item['name']).read_bytes()
            if sha(raw)!=item['sha256'] or len(raw)!=item['size']:raise ValueError('changed evidence part')
            digest.update(raw);size+=len(raw);stream.write(raw)
    if size!=manifest['archive_size'] or digest.hexdigest()!=manifest['archive_sha256']:
        raise ValueError('changed joined evidence')
    evidence=output/'evidence';extract_checked(archive,evidence,manifest['entries'])
    sources=json.loads((evidence/'source-snapshot.json').read_text(encoding='utf-8'))
    with zipfile.ZipFile(evidence/'source-snapshot.zip') as z:
        sizes={i.filename:i.file_size for i in z.infolist()}
    extract_checked(evidence/'source-snapshot.zip',evidence/'source',
                    {name:dict(sha256=digest,size=sizes[name]) for name,digest in sources.items()})
    script=evidence/'source/interim'/HERE.name/'measure.py'
    measured=output/'material-history.json'
    subprocess.run([sys.executable,'-I','-X','utf8',str(script),'--data',str(evidence/'data'),'--output',str(measured)],check=True)
    if measured.read_bytes()!=(evidence/'expected-material-history.json').read_bytes():
        raise ValueError('sealed remeasurement differs; no acceptance or relaxed comparison')
    # Whole saved Run equality, including states, input identity, continuation,
    # all trial reports, and all frames. Only measured solver time is excluded.
    code="""import json,sys
from pathlib import Path
root=Path(sys.argv[1]);sys.path.insert(0,str(root/'source/tools'))
from probe_low_preload_nc import deterministic
read=lambda p: json.loads(p.read_text(encoding='utf-8'))
data=root/'data'
for n in (128,256):
 direct=deterministic(read(data/f'direct{n}/run.json'))
 assert deterministic(read(data/f'trace{n}/run.json'))==direct
 for stage in (2,3):
  assert deterministic(read(data/f'restarts/n{n}-stage{stage}/run.json'))==direct
"""
    subprocess.run([sys.executable,'-I','-X','utf8','-c',code,str(evidence)],check=True)
    report=dict(all_parts_and_entries_verified=True,parts=len(manifest['parts']),entries=len(manifest['entries']),
        core_source_files_verified=len(sources),sealed_sources_isolated_utf8=True,
        material_history_byte_identical=True,whole_saved_runs_exact_except_elapsed_seconds=True,
        numerical_acceptance=False,completion_restart_verified=False)
    (output/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2));return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
