"""Bind geostatic/slope excavation evidence to actual inputs, sources and binary."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def record(run,rust_log,lint_log,format_log):
    root=Path(__file__).resolve().parents[1];solver=run/'solver/soil-fem-3d.exe';digest=sha(solver)
    suite=json.loads((run/'verification.json').read_text());assert suite['passed'] and suite['solver_sha256']==digest
    assert sha(root/'target/release/soil-fem-3d.exe')==digest
    text=rust_log.read_text();count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',text)))
    assert count==69
    for p in (rust_log,lint_log,format_log):assert p.read_text().rstrip().endswith('Exit code: 0')
    cad=(run/'cad-tests.log').read_text();assert 'Ran 5 tests' in cad and 'OK' in cad and cad.rstrip().endswith('Exit code: 0')
    independent=json.loads((run/'independent-slope.json').read_text());assert independent['solver_sha256']==digest
    assert len(independent['frames'])==3 and max(f['relative_3d_displacement_error'] for f in independent['frames'])<1e-9
    pit_reference=json.loads((run/'independent-pit.json').read_text());assert pit_reference['solver_sha256']==digest
    assert len(pit_reference['frames'])==2 and max(f['relative_3d_displacement_error'] for f in pit_reference['frames'])<1e-9
    assert pit_reference['frames'][1]['active_elements']<pit_reference['frames'][0]['active_elements']
    metrics={};files=[solver,rust_log,lint_log,format_log,run/'cad-tests.log',run/'verification.json',run/'independent-slope.json',run/'independent-pit.json']
    for name in ('mc-slope','mcc-slope','free-mcc','slope-pit/analysis','mc-pit','mcc-regression-cad','initial-only'):
        folder=run/name;package_path=folder/'accepted-stage-run.json';package=json.loads(package_path.read_text())
        assert package['solver_binary_sha256']==sha(folder/'solver/soil-fem-3d.exe')==digest
        for path,expected in package['inputs_and_results_sha256'].items():
            p=Path(path);assert sha(p)==expected;files.append(p)
        result_path=folder/('initial-equilibrium.json' if name=='initial-only' else 'stage-results.json')
        result=json.loads(result_path.read_text());assert result['complete']
        summary={'frames':len(result['frames']),'integration_points':len(result['integration_points'])}
        if result.get('geostatic'):
            g=result['geostatic'];assert g['complete'] and g['final_relative_physical_force_residual']<1e-9
            summary.update(initial_force_residual=g['final_relative_physical_force_residual'],
                initial_maximum_strain=g['maximum_absolute_engineering_strain_from_seed'],
                initial_plastic_points=sum(any(abs(x)>1e-9 for x in p['plastic_strain_compression_positive']) for p in g['last_committed']['material_states']))
        if name=='slope-pit/analysis':
            mesh=json.loads((folder/'ground/mesh.json').read_text());xyz=mesh['nodes']
            surface=mesh['node_sets']['Pit:surrounding_surface']
            assert surface and min(xyz[n][2] for n in surface)<max(p[2] for p in xyz)-1.
            assert len(result['excavation_audits'])==1
            audit=result['excavation_audits'][0];assert audit['zero_release_relative_force_residual']<1e-8
            assert any(abs(x)>1e-9 for x in audit['initial_holding_load_kn'])
            assert all(x==0 for x in result['last_committed']['excavation_release_load_kn'])
            rejected=[a for a in result['attempts'] if not a['accepted']]
            assert rejected and all(a['committed_state_unchanged_on_rejection'] for a in rejected)
            summary['final_monitors']=result['frames'][-1]['monitors']
            summary['active_elements']=sum(result['last_committed']['active_elements'])
            summary['removed_elements']=len(result['last_committed']['active_elements'])-summary['active_elements']
        metrics[name]=summary;files.extend([package_path,folder/'stage-verification.json'])
    for name in ('checkpoint','slope-pit-checkpoint'):
        p=run/name/'verification.json';checkpoint=json.loads(p.read_text());assert checkpoint['solver_sha256']==digest
        assert len(checkpoint['restarts'])==3 and all(r['exact_histories_frames_reference_and_increment_path'] for r in checkpoint['restarts'])
        assert checkpoint['modified_holding_load_rejected'] and checkpoint['strain_limit_preserves_seed_and_rejects_export']
        for path,expected in checkpoint['artifacts_sha256'].items():assert sha(run/name/path)==expected
        files.append(p)
    initial=run/'initial-only/accepted-run.json';accepted=json.loads(initial.read_text());assert accepted['equilibrium_verified']
    for name,expected in accepted['files_sha256'].items():assert sha(initial.parent/name)==expected
    regression=run/'mcc-regression.json';compatible=json.loads(regression.read_text());assert compatible['exact_committed_state'] and compatible['exact_frames']
    files.extend([initial,regression,run/'free-mcc-input.seed-provenance.json',run/'slope-pit/ground.seed-provenance.json'])
    files.extend(root/name for name in ('Cargo.lock','Cargo.toml','sources.lock.json','requirements.lock.txt','rust-toolchain.toml','README.md'))
    files.extend(p for name in ('src','tools','scripts','tests','examples','docs') for p in (root/name).rglob('*')
        if p.is_file() and p.suffix in ('.rs','.py','.ps1','.json','.md'))
    report={'schema':'soil-fem-geostatic-final-validation/1','passed':True,'rust_tests':count,'cad_tests':5,
        'solver_sha256':digest,'metrics':metrics,'multi_phase_restart_exact':True,'mcc_legacy_exact':True,
        'independent_reference_scope':'Full 3D elastic MC slope and excavation solves; nonlinear MCC force/moment/state audits and prior constitutive/column references, not an independent nonlinear 3D plastic displacement solve',
        'scope':'General initial state and slope/local excavation milestone; P6 overall incomplete',
        'sha256':{str(p.relative_to(root)):sha(p) for p in sorted(set(files))}}
    with (run/'final-validation.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--rust-log',type=Path,required=True)
    p.add_argument('--lint-log',type=Path,required=True);p.add_argument('--format-log',type=Path,required=True)
    a=p.parse_args();record(*[v.resolve() for v in (a.run,a.rust_log,a.lint_log,a.format_log)])
