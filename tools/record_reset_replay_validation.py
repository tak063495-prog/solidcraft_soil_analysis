"""Bind RESET, fresh replay, rejected runs, independent solves and old results."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import meshio
import numpy as np


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def record(run,rust,lint,format_log):
    root=Path(__file__).resolve().parents[1];solver=run/'solver/soil-fem-3d.exe';digest=sha(solver)
    suite=json.loads((run/'verification.json').read_text());assert suite['passed'] and suite['solver_sha256']==digest==sha(root/'target/release/soil-fem-3d.exe')
    count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust.read_text())));assert count==82
    for p in (rust,lint,format_log):assert p.read_text().rstrip().endswith('Exit code: 0')
    reference=json.loads((run/'independent.json').read_text());assert reference['solver_sha256']==digest and len(reference['frames'])==11
    assert max(f['relative_3d_displacement_error'] for f in reference['frames'])<1e-9
    checkpoint=json.loads((run/'checkpoint/verification.json').read_text());assert checkpoint['solver_sha256']==digest and len(checkpoint['restarts'])==11
    assert all(c['exact_states_frames_material_ids_and_audits'] for c in checkpoint['restarts'])
    files=[solver,rust,lint,format_log,run/'verification.json',run/'independent.json',run/'checkpoint/verification.json']
    for name,expected in checkpoint['sha256'].items():p=run/'checkpoint'/name;assert sha(p)==expected;files.append(p)
    replay=json.loads((run/'replay-verification.json').read_text());assert replay['solver_sha256']==digest and replay['unit_factor_exact'] and replay['stability_classification']=='not_evaluated'
    assert len(replay['trials'])==3 and max(t['maximum_independent_3d_displacement_error'] for t in replay['trials'])<1e-9
    for name,expected in replay['sha256'].items():p=Path(name);assert sha(p)==expected;files.append(p)
    files.append(run/'replay-verification.json')
    failures=run/('failure-guards' if (run/'failure-guards').exists() else 'failures');failure=json.loads((failures/'verification.json').read_text())
    assert failure['solver_sha256']==digest and len(failure['cases'])==2
    assert all(not c['noncompleted_accepted_exports'] and not c['stability_or_fs_inferred'] for c in failure['cases'])
    for name,expected in failure['sha256'].items():p=failures/name;assert sha(p)==expected;files.append(p)
    files.append(failures/'verification.json')
    for name in ('numerical','policy'):assert (failures/name/'replay.log').read_text().rstrip().endswith('Exit code: 1');files.append(failures/name/'replay.log')
    regression=run/'geostatic-regression';geo=json.loads((regression/'verification.json').read_text());assert geo['passed'] and geo['solver_sha256']==digest
    cad=regression/'cad-tests.log';assert 'Ran 5 tests' in cad.read_text() and cad.read_text().rstrip().endswith('Exit code: 0')
    compat=regression/'mcc-regression.json';compatibility=json.loads(compat.read_text());assert compatibility['exact_committed_state'] and compatibility['exact_frames']
    files.extend([regression/'verification.json',cad,compat])
    folders=[run/'cad']+[run/f'replay/trials/trial-{i:03}' for i in range(3)]+[regression/name for name in ('mc-slope','mcc-slope','free-mcc','slope-pit/analysis','mc-pit','mcc-regression-cad','initial-only')]
    for folder in folders:
        p=folder/'accepted-stage-run.json';accepted=json.loads(p.read_text());assert accepted['solver_binary_sha256']==digest
        for path,expected in accepted['inputs_and_results_sha256'].items():f=Path(path);assert sha(f)==expected;files.append(f)
        files.extend([p,folder/'stage-verification.json'])
    data=json.loads((run/'cad/stage-results.json').read_text());mesh=json.loads((run/'cad/ground/mesh.json').read_text());assert data['complete'] and data['schema']=='soil-fem-stage-result/6'
    for i,frame in enumerate(data['frames']):
        state=frame['state'];vtk=run/f'cad/results/stage-{i:03}.vtu';v=meshio.read(vtk);used=v.point_data['original_node_id'].astype(int)
        u=np.asarray(state['displacement_m']);display_reference=np.asarray(state.get('display_reference_displacement_m') or data['reference_displacement_m'])
        assert np.array_equal(v.point_data['display_reference_displacement_m'],display_reference[used]);assert np.allclose(v.point_data['displacement_m'],(u-display_reference)[used],atol=1e-14,rtol=1e-12)
        if i and data['case']['stages'][i-1]['kind']=='reset_display':
            before=data['frames'][i-1]['state'];assert state['material_states']==before['material_states'] and state['displacement_m']==before['displacement_m'] and state['applied_nodal_load_kn']==before['applied_nodal_load_kn']
        files.append(vtk)
    files.extend(root/p for p in ('Cargo.lock','Cargo.toml','requirements.lock.txt','sources.lock.json','rust-toolchain.toml','README.md'))
    files.extend(p for name in ('src','tools','scripts','tests','examples','docs') for p in (root/name).rglob('*') if p.is_file() and p.suffix in ('.rs','.py','.ps1','.json','.md'))
    report={'schema':'soil-fem-reset-replay-final-validation/1','passed':True,'rust_tests':count,'cad_tests':5,'solver_sha256':digest,
        'nodes':len(mesh['nodes']),'tetrahedra10':len(mesh['tetrahedra10']),'integration_points':len(data['integration_points']),
        'separate_process_exact_restarts':11,'independent_replayed_factors':[t['factor'] for t in replay['trials']],
        'maximum_independent_3d_displacement_error':max(f['relative_3d_displacement_error'] for f in reference['frames']),
        'maximum_absolute_display_displacement_error_m':max(f['absolute_display_displacement_error_m'] for f in reference['frames']),
        'unit_factor_exact':True,'mcc_legacy_exact':True,'actual_cli_failure_and_policy_rejection_retained':True,
        'remaining':'2D deformed-reference RESET compatibility, selected intermediate replay origins, physical Fs search/gates and P6 overall remain incomplete',
        'sha256':{str(p.relative_to(root)):sha(p) for p in sorted(set(files))}}
    with (run/'final-validation.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--rust-log',type=Path,required=True);p.add_argument('--lint-log',type=Path,required=True);p.add_argument('--format-log',type=Path,required=True)
    a=p.parse_args();record(*[v.resolve() for v in (a.run,a.rust_log,a.lint_log,a.format_log)])
