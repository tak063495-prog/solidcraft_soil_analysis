"""Bind jacking/slack, independent soil solve, recovery and regression evidence."""
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
    if (run/'zero-lock-guard.json').exists():raise ValueError('superseded zero-lock run cannot be accepted')
    solver=run/'solver/soil-fem-3d.exe';digest=sha(solver);suite=read(run/'verification.json')
    assert suite['passed'] and suite['solver_sha256']==digest
    for p in checks.values():assert p.read_text().rstrip().endswith('Exit code: 0')
    rust=checks['rust'].read_text();count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)));assert count>=98
    for name in ('curved_tet10_inverse_and_interface_hosts_keep_coordinate_and_activity',
                 'anchor_jack_lock_slack_and_plastic_soil_restore_exactly_on_cutback_and_restart',
                 'excavation_with_retained_anchor_and_beam_releases_only_removed_soil_stress'):
        assert f'test {name} ... ok' in rust
    cad=(run/'cad-tests.log').read_text();assert 'Ran 5 tests' in cad and cad.rstrip().endswith('Exit code: 0')
    examples={};files=list(checks.values());receipts=[];restarts=[]
    for name in ('cad','mcc'):
        folder=run/name;data=read(folder/'stage-results.json');mesh=read(folder/'ground/mesh.json')
        assert data['complete'] and data['schema']=='soil-fem-stage-result/8'
        accepted=read(folder/'accepted-stage-run.json');assert accepted['solver_binary_sha256']==digest==sha(folder/'solver/soil-fem-3d.exe')
        receipts.append(accepted['inputs_and_results_sha256'])
        active=[len(f['anchor_results']['anchors']) for f in data['frames']];assert active==[0,1,1,1,1,1,0]
        forces=[]
        for fi,f in enumerate(data['frames']):
            anchors=f['anchor_results']['anchors'];forces.append(anchors[0]['segment_forces_kn'][0] if anchors else None)
            p=folder/f'results/anchors-{fi:03}.vtu'
            if anchors:
                a=anchors[0];assert abs(a['free_length_m']-5)<1e-10 and abs(a['bonded_length_m']-5)<1e-10
                vtk=meshio.read(p);assert len(vtk.cells_dict['line'])==4
                assert np.allclose(vtk.cell_data['tension_kn'][0],a['segment_forces_kn'])
                assert np.allclose(vtk.point_data['total_displacement_from_cad_m'],a['point_displacement_m'])
                assert np.allclose(vtk.point_data['force_from_anchor_to_soil_kn'],a['point_soil_force_kn'])
            else:
                p=ET.parse(p).find('UnstructuredGrid/Piece');assert p.get('NumberOfCells')=='0'
        assert abs(forces[2]-50)<1e-8 and forces[3]==0 and abs(forces[5])<1e-8
        if name=='cad':assert abs(forces[4]-50)<1e-8
        else:
            assert forces[4]==0
            before=data['frames'][4]['state']['anchors']['anchors'][0]
            after=data['frames'][5]['state']['anchors']['anchors'][0]
            assert before['reference_extensions_m'][0]==after['reference_extensions_m'][0]
            assert not before['segment_taut'][0] and not after['segment_taut'][0]
        verify=read(folder/'stage-verification.json')
        assert max(f['independent_force_residual'] for f in verify['frames'])<1e-8
        assert max(f['relative_moment_error'] for f in verify['frames'])<1e-7
        h=read(run/('mcc-history' if name=='mcc' else 'history')/'verification.json')
        assert h['solver_sha256']==digest and h['exact_final_and_all_frames'] and h['rejected_structure_and_material_changed']
        if name=='mcc':assert h['rejected_plastic_changed']
        receipts.append(h['sha256'])
        cpname='mcc-checkpoint' if name=='mcc' else 'checkpoint';cp=read(run/cpname/'verification.json')
        assert cp['solver_sha256']==digest and len(cp['restarts'])==7 and any(r['saved_anchor_force_control'] for r in cp['restarts'])
        assert all(r['exact_states_frames_material_ids_and_audits'] for r in cp['restarts']);restarts+=cp['restarts']
        for path,expected in cp['sha256'].items():assert sha(run/cpname/path)==expected
        examples[name]={'nodes':len(mesh['nodes']),'tetrahedra10':len(mesh['tetrahedra10']),'integration_points':len(data['integration_points']),
                        'active_anchor_counts':active,'head_tension_kn_by_frame':forces,
                        'maximum_independent_force_residual':max(f['independent_force_residual'] for f in verify['frames'])}
    reference=read(run/'independent.json');assert reference['solver_sha256']==digest and len(reference['frames'])==7
    assert max(f['relative_3d_displacement_error'] for f in reference['frames'])<1e-7
    receipts.append(reference['sha256'])
    for receipt in receipts:
        for path,expected in receipt.items():assert sha(Path(path))==expected,path
    regression=run/'beam-regression';b=read(regression/'verification.json');assert b['passed'] and b['solver_sha256']==digest
    bcp=read(regression/'mcc-checkpoint/verification.json');assert len(bcp['restarts'])==5 and all(r['exact_states_frames_material_ids_and_audits'] for r in bcp['restarts'])
    geo=regression/'geostatic-regression';g=read(geo/'verification.json');assert g['passed'] and g['solver_sha256']==digest
    compat=geo/'mcc-regression.json'
    if compat.exists():
        c=read(compat);assert c['exact_committed_state'] and c['exact_frames']
        for path,expected in c['files_sha256'].items():p=Path(path);assert sha(p)==expected;files.append(p)
    for name in ('cad','independent','history','checkpoint','mcc','mcc-history','mcc-checkpoint','beam-regression'):
        assert (run/(name+'.log')).read_text().rstrip().endswith('Exit code: 0')
    for path in run.rglob('accepted-stage-run.json'):
        package=read(path);assert package['solver_binary_sha256']==digest
        for name,expected in package['inputs_and_results_sha256'].items():assert sha(Path(name))==expected
    for path in (regression/'independent.json',regression/'existing-independent.json'):
        comparison=read(path);assert comparison['solver_sha256']==digest
        assert max(f['relative_3d_displacement_error'] for f in comparison['frames'])<1e-7
        for name,expected in comparison['sha256'].items():assert sha(Path(name))==expected
    sources=[root/name for name in ('Cargo.lock','Cargo.toml','rust-toolchain.toml','requirements.lock.txt','sources.lock.json','README.md')]
    sources += [p for name in ('src','tests','scripts','examples','docs') for p in (root/name).rglob('*') if p.is_file() and p.suffix in ('.rs','.py','.ps1','.json','.md')]
    sources += list((root/'tools').glob('*.py'));source_map={p.relative_to(root).as_posix():sha(p) for p in sorted(set(sources))}
    with zipfile.ZipFile(run/'source-snapshot.zip','x',compression=zipfile.ZIP_DEFLATED) as archive:
        for name,expected in source_map.items():
            payload=(root/name).read_bytes();assert hashlib.sha256(payload).hexdigest()==expected;archive.writestr(name,payload)
    (run/'source-snapshot.json').write_text(json.dumps(source_map,indent=2));files += [p for p in run.rglob('*') if p.is_file()]
    report={'schema':'soil-fem-anchors-final-validation/1','passed':True,'solver_sha256':digest,'rust_tests':count,'cad_tests':5,'examples':examples,
        'separate_process_exact_restarts':len(restarts),'jacking_force_control_checkpoint_verified':True,
        'maximum_independent_3d_displacement_error':max(f['relative_3d_displacement_error'] for f in reference['frames']),
        'maximum_independent_tension_error_kn':max(f['maximum_tension_error_kn'] for f in reference['frames']),
        'mc_and_mcc_cutback_exact_vs_clean':True,'legacy_mcc_exact':compat.exists(),
        'zero_force_relock_preserves_existing_slack':True,
        'limitations':'Straight massless elastic tension-only tendons; perfectly bonded segmented anchorage, no bond-slip/pullout capacity. Head attaches to a continuum host. Global independent reference covers elastic MC soil; plastic MCC is audited and replayed, not independently certified. P6 overall incomplete.',
        'sha256':{p.relative_to(root).as_posix():sha(p) for p in sorted(set(files))}}
    with output.open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True)
    for name in ('rust','format','clippy','build'):p.add_argument('--'+name+'-log',type=Path)
    a=p.parse_args();run=a.run.resolve();defaults={'rust':'rust-tests.log','format':'format.log','clippy':'lint.log','build':'build.log'}
    record(run,{name:(getattr(a,name+'_log') or run/path).resolve() for name,path in defaults.items()})
