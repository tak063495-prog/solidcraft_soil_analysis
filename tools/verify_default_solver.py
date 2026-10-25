"""Bind the normal target/release build to accepted mixed and MCC results.

Different Cargo target directories can produce different executable hashes.
Do not substitute one executable for another to force a hash match. Execute
the actual default binary and compare every result field, excluding only the
declared executable fingerprint and measured wall-clock solve durations.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from freeze_sources import freeze
from verify_acceptance_packet import sha,verify as verify_packet


def read(path):return json.loads(path.read_text(encoding='utf-8'))


def deterministic(value):
    if isinstance(value,dict):
        return {k:deterministic(v) for k,v in value.items() if k not in ('solver_binary_sha256','elapsed_seconds')}
    if isinstance(value,list):return [deterministic(v) for v in value]
    return value


def run(out):
    root=Path(__file__).resolve().parents[1];solver=root/'target/release/soil-fem-3d.exe'
    accepted=root/'runs/mixed-mechanics-acceptance-20261025-v1'
    cycle=root/'runs/mcc-tet10-cycle-acceptance-20261025-v1'
    parent_integrity=verify_packet(accepted);cycle_integrity=verify_packet(cycle)
    binary_hash=sha(solver);sources=read(accepted/'source-snapshot.json')
    for name,digest in sources.items():
        if name.startswith('src/') or name in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):
            assert sha(root/name)==digest,name
    build=root/'runs/mixed-mechanics-production-build-20261025-v1.log'
    assert build.read_text().rstrip().endswith('Exit code: 0')
    out.mkdir(exist_ok=False);(out/'solver').mkdir();shutil.copyfile(solver,out/'solver/soil-fem-3d.exe')
    cases=[]
    for n in (2,4,6):
        for nu in ('0_3','0_49','0_4999'):
            tag=f'tet10-n{n}-nu{nu}';base=accepted/'mixed-data'
            cases.append((tag,'solve-mixed-mechanics',[base/f'tet10-n{n}-mesh.json',base/(tag+'-case.json')],base/(tag+'-result.json')))
    for name in ('column','rotated'):
        base=accepted/'column-data';cases.append((name,'solve-mixed-mechanics',[base/(name+'-mesh.json'),base/(name+'-case.json')],base/(name+'-result.json')))
    base=accepted/'curved-cad-data'
    cases.append(('curved-cad-q5','solve-mixed-mechanics',[base/'mesh-q5.json',base/'case.json'],base/'result-q5.json'))
    cases.append(('mcc-load-unload-cut','solve-stages',[cycle/'inputs/mesh.json',cycle/'inputs/ground.json',cycle/'inputs/slope-cycle-case.json'],accepted/'current-mcc/result.json'))
    inputs=[solver,Path(__file__),root/'tools/run_logged.py',build]
    inputs.extend(p for _,_,paths,expected in cases for p in [*paths,expected])
    fingerprints={str(p):sha(p) for p in set(inputs)};reports=[]
    for name,command,paths,expected_path in cases:
        result=out/(name+'-result.json')
        subprocess.run([sys.executable,str(root/'tools/run_logged.py'),str(out/(name+'.log')),str(solver),command,
                        *map(str,paths),str(result)],cwd=root,check=True)
        actual=read(result);expected=read(expected_path)
        fingerprint=actual.get('solver_binary_sha256',(actual.get('provenance') or {}).get('solver_binary_sha256'))
        assert fingerprint==binary_hash,name
        left=json.dumps(deterministic(actual),sort_keys=True,separators=(',',':'),allow_nan=False)
        right=json.dumps(deterministic(expected),sort_keys=True,separators=(',',':'),allow_nan=False)
        assert left==right,'non-runtime result differs: '+name
        reports.append(dict(name=name,all_nonruntime_result_fields_exact=True,reference_sha256=sha(expected_path),actual_sha256=sha(result)))
    assert fingerprints=={str(p):sha(p) for p in set(inputs)},'default executable or input changed during verification'
    freeze(out)
    record=dict(schema='soil-fem-default-solver-equivalence/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Actual default build: 9 mixed cubes, two-layer/rotated columns, curved CAD, and legacy MCC cycle; exact accepted-field comparison',
        solver_sha256=binary_hash,accepted_mixed_solver_sha256=parent_integrity['solver_sha256'],
        parent_manifest_sha256=parent_integrity['manifest_sha256'],mcc_parent_manifest_sha256=cycle_integrity['manifest_sha256'],
        cases=reports,excluded_keys=['solver_binary_sha256','elapsed_seconds'],
        source_snapshot_sha256=sha(out/'source-snapshot.zip'),input_and_source_sha256=fingerprints,
        limitations='Build equivalence for existing accepted cases, not additional model or full P6 certification',
        sha256={p.relative_to(root).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'final-validation.json').open('x') as f:json.dump(record,f,indent=2)
    print(json.dumps(dict(passed=True,cases=len(reports),all_nonruntime_result_fields_exact=True,solver_sha256=binary_hash)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();run(a.output.resolve())
