"""Reuse successful Rust checks only for byte-identical compilation/test inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import zipfile
from freeze_sources import manifest


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))


def reuse(previous,current):
    before=read(previous/'source-snapshot.json');after=read(current/'source-snapshot.json')
    assert after==manifest(),'current source changed after freeze'
    def compilation_input(name):
        return name in ('Cargo.toml','Cargo.lock','rust-toolchain.toml') or name.startswith(('src/','examples/')) or (name.startswith('tests/') and name.endswith('.rs'))
    a={n:h for n,h in before.items() if compilation_input(n)}
    b={n:h for n,h in after.items() if compilation_input(n)}
    assert a==b and a,'Rust source, dependency, example or test input changed'
    with zipfile.ZipFile(previous/'source-snapshot.zip') as archive:
        assert archive.testzip() is None
        for name,expected in before.items():assert hashlib.sha256(archive.read(name)).hexdigest()==expected
    logs={}
    for name in ('rust-tests','format','lint','build'):
        path=previous/(name+'.log');content=path.read_text(encoding='utf-8')
        assert content.rstrip().endswith('Exit code: 0'),name
        logs[str(path.resolve())]=sha(path)
    rust=(previous/'rust-tests.log').read_text(encoding='utf-8')
    assert sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)))>=121
    shutil.copyfile(previous/'rust-tests.log',current/'rust-tests.log')
    proof=current/'reused-rust-check';proof.mkdir()
    for name in ('source-snapshot.json','source-snapshot.zip'):
        shutil.copyfile(previous/name,proof/name)
    report={'schema':'soil-fem-reused-rust-test-evidence/1','passed':True,
        'previous_run':str(previous),'expected_solver_sha256':sha(previous/'solver/soil-fem-3d.exe'),
        'identical_compilation_and_test_inputs':a,'original_check_logs_sha256':logs,
        'previous_source_snapshot_sha256':sha(previous/'source-snapshot.zip'),
        'scope':'Rust test evidence reused from a run whose later integration verification failed. All Rust compilation and test inputs match exactly; format, lint, build and integrations are run again.'}
    with (current/'rust-test-reuse.json').open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print('Reused successful Rust test evidence for',len(a),'identical compilation/test input files')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('previous',type=Path);p.add_argument('current',type=Path)
    a=p.parse_args();reuse(a.previous.resolve(),a.current.resolve())
