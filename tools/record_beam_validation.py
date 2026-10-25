"""Freeze beam integration evidence without promoting the overall P6 goal."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile
import meshio
import numpy as np


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))


def record(run,checks):
    root=Path(__file__).resolve().parents[1];output=run/'final-validation.json'
    if output.exists():raise ValueError('validation exists')
    solver=run/'solver/soil-fem-3d.exe';digest=sha(solver);suite=read(run/'verification.json')
    assert suite['passed'] and suite['solver_sha256']==digest
    for p in checks.values():assert p.read_text().rstrip().endswith('Exit code: 0')
    rust=checks['rust'].read_text();count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)));assert count>=89
    for name in ('timoshenko_uniform_load_and_nonzero_constraint_lifting_match_theory',
                 'rehashed_checkpoint_rejects_corrupted_recovered_beam_kinematics',
                 'condensed_beam_is_coupled_to_mc_and_mcc_force_moment_and_exact_restart'):
        assert f'test {name} ... ok' in rust
    cad=(run/'cad-tests.log').read_text();assert 'Ran 5 tests' in cad and cad.rstrip().endswith('Exit code: 0')
    examples={};files=list(checks.values());source_receipts=[]
    for name in ('cad','mcc','existing'):
        folder=run/name;data=read(folder/'stage-results.json');mesh=read(folder/'ground/mesh.json')
        assert data['complete'] and data['schema']=='soil-fem-stage-result/7'
        accepted=read(folder/'accepted-stage-run.json');assert accepted['solver_binary_sha256']==digest
        for path,expected in accepted['inputs_and_results_sha256'].items():assert sha(Path(path))==expected
        assert sha(folder/'solver/soil-fem-3d.exe')==digest
        active=[f['beam_results']['active_members'] for f in data['frames']]
        assert active==([2,2,2,0] if name=='existing' else [0,2,2,2,0])
        for fi,frame in enumerate(data['frames']):
            b=frame['beam_results'];assert b['true_relative_interior_residual']<1e-9
            assert max(map(abs,frame['global_moment_balance_error_kn_m']))<1e-6
            p=folder/f'results/beams-{fi:03}.vtu'
            if b['active_members']:
                vtk=meshio.read(p);assert len(vtk.cells_dict['line'])==b['active_members']
                used=vtk.point_data['node_index'].astype(int)
                assert np.allclose(vtk.point_data['total_displacement_from_cad_m'],np.asarray(b['node_displacement_m'])[used])
                assert np.allclose(vtk.point_data['total_rotation_rad'],np.asarray(b['node_rotation_rad'])[used])
            else:
                piece=ET.parse(p).find('UnstructuredGrid/Piece');assert piece.get('NumberOfCells')==piece.get('NumberOfPoints')=='0'
        if name!='existing':
            h=read(run/('mcc-history' if name=='mcc' else 'history')/'verification.json')
            assert h['exact_final_and_all_frames'] and h['rejected_structure_and_material_changed'] and h['solver_sha256']==digest
            if name=='mcc':assert h['rejected_plastic_changed']
            source_receipts.append(h['sha256'])
        examples[name]={'nodes':len(mesh['nodes']),'tetrahedra10':len(mesh['tetrahedra10']),
                        'integration_points':len(data['integration_points']),
                        'active_member_counts':active,'maximum_beam_interior_residual':max(f['beam_results']['true_relative_interior_residual'] for f in data['frames'])}
    comparisons=[read(run/'independent.json'),read(run/'existing-independent.json')]
    for c in comparisons:
        assert c['solver_sha256']==digest
        assert max(f['relative_3d_displacement_error'] for f in c['frames'])<1e-7
        source_receipts.append(c['sha256'])
    restarts=[]
    for name in ('checkpoint','mcc-checkpoint'):
        c=read(run/name/'verification.json');assert c['solver_sha256']==digest and len(c['restarts'])==5
        assert all(r['exact_states_frames_material_ids_and_audits'] for r in c['restarts']);restarts+=c['restarts']
        for path,expected in c['sha256'].items():assert sha(run/name/path)==expected
    for receipt in source_receipts:
        for path,expected in receipt.items():assert sha(Path(path))==expected
    regression=run/'geostatic-regression';g=read(regression/'verification.json');assert g['passed'] and g['solver_sha256']==digest
    compat=regression/'mcc-regression.json'
    if compat.exists():
        r=read(compat);assert r['exact_committed_state'] and r['exact_frames']
        for path,expected in r['files_sha256'].items():
            p=Path(path);assert sha(p)==expected;files.append(p)
    for name in ('cad','independent','history','checkpoint','mcc','mcc-history','mcc-checkpoint','existing','existing-independent','geostatic-regression'):
        assert (run/(name+'.log')).read_text().rstrip().endswith('Exit code: 0')
    sources=[root/name for name in ('Cargo.lock','Cargo.toml','rust-toolchain.toml','requirements.lock.txt','sources.lock.json','README.md')]
    sources += [p for name in ('src','tests','scripts','examples','docs') for p in (root/name).rglob('*') if p.is_file() and p.suffix in ('.rs','.py','.ps1','.json','.md')]
    sources += list((root/'tools').glob('*.py'));source_map={p.relative_to(root).as_posix():sha(p) for p in sorted(set(sources))}
    with zipfile.ZipFile(run/'source-snapshot.zip','x',compression=zipfile.ZIP_DEFLATED) as archive:
        for name,expected in source_map.items():
            payload=(root/name).read_bytes();assert hashlib.sha256(payload).hexdigest()==expected;archive.writestr(name,payload)
    (run/'source-snapshot.json').write_text(json.dumps(source_map,indent=2))
    files += [p for p in run.rglob('*') if p.is_file()]
    report={'schema':'soil-fem-beams-final-validation/1','passed':True,'solver_sha256':digest,'rust_tests':count,'cad_tests':5,
            'examples':examples,'separate_process_exact_restarts':len(restarts),
            'mc_and_mcc_cutback_exact_vs_clean_accepted_path':True,
            'maximum_independent_3d_displacement_error':max(f['relative_3d_displacement_error'] for c in comparisons for f in c['frames']),
            'maximum_independent_end_action_error':max(f['maximum_end_action_error'] for c in comparisons for f in c['frames']),
            'legacy_mcc_exact':compat.exists(),
            'limitations':'Linear small-strain frames; soil attaches by shared translations, rotations have explicit external supports. No soil moment bond, embedded interpolation, beam plasticity/buckling, anchors, plates, history JOINT or P5/P6 completion.',
            'sha256':{p.relative_to(root).as_posix():sha(p) for p in sorted(set(files))}}
    with output.open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True)
    for name in ('rust','format','clippy','build'):p.add_argument('--'+name+'-log',type=Path)
    a=p.parse_args();run=a.run.resolve();defaults={'rust':'rust-tests.log','format':'format.log','clippy':'lint.log','build':'build.log'}
    record(run,{name:(getattr(a,name+'_log') or run/path).resolve() for name,path in defaults.items()})
