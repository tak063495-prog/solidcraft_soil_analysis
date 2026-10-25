"""Seal fixed-mesh projected-MCC integration evidence, not general/P6 acceptance."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import zipfile
from freeze_sources import freeze


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def receipt(path):
    value=read(path);assert value['passed'],path
    for name,digest in value.get('input_and_source_sha256',{}).items():assert sha(Path(name))==digest,name
    for name,digest in value.get('evidence_sha256',{}).items():assert sha(path.parent/name)==digest,name
    return value


def record(study,assessment_dir,output):
    root=Path(__file__).resolve().parents[1];runs=root/'runs'
    assessment=read(assessment_dir/'assessment.json')
    assert assessment['integration_gates_passed'] and not assessment['production_backend_accepted']
    assert assessment['reference_refinement_includes_initialization'] and not assessment['spatial_convergence_verified']
    assert read(study/'reference-400-progress.json')['state']=='completed'
    assert (runs/'projected-mcc-reference-400-driver-20261025-v1.log').read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
    for name,digest in assessment['sha256'].items():assert sha(Path(name))==digest,name
    lifecycle=runs/'projected-mcc-lifecycle-20261025-v3'
    exported=runs/'projected-mcc-export-verification-20261025-v1'
    legacy_export=runs/'projected-mcc-legacy-export-regression-20261025-v1'
    legacy_solver=runs/'projected-mcc-legacy-solver-regression-20261025-v1'
    life,ex,le,ls=[receipt(p/'verification.json') for p in (lifecycle,exported,legacy_export,legacy_solver)]
    assert len(life['restarts'])==3 and life['rejected_converged_trials']==3 and life['rejected_plastic_and_domain_changes']
    assert life['full_control_frames_and_states_exact'] and life['changed_load_checkpoint_refused']
    assert life['restarts'][-1]['maximum_frozen_compatibility_offset_m3']>1e-5
    assert len(ex['rejected_inputs'])==14 and len(le['rejected_inputs'])==15
    assert len(ls['cases'])==2 and all(c['all_nonruntime_fields_exact'] for c in ls['cases'])
    binary=study/'solver/soil-fem-3d.exe';digest=sha(binary)
    assert digest=='09eed93197f8aedcbae79b9d404ccab6becb8b5a0d8b95135cef623a9578210c'
    assert ls['solver_sha256']==digest
    logs=['projected-mcc-release-tests-20261025-v1.log','projected-mcc-clippy-20261025-v1.log',
          'projected-mcc-fmt-20261025-v1.log','projected-mcc-build-20261025-v1.log',
          'projected-mcc-lifecycle-20261025-v3.log','projected-mcc-export-verification-20261025-v1.log',
          'projected-mcc-legacy-export-regression-20261025-v1.log','projected-mcc-legacy-solver-regression-20261025-v1.log']
    for name in logs:assert (runs/name).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    counts=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored',(runs/logs[0]).read_text())
    assert sum(int(n) for n,_,_ in counts)==242 and all(f==i=='0' for _,f,i in counts)
    build=study/'build-sources';sources=read(build/'source-snapshot.json');assert len(sources)==488
    with zipfile.ZipFile(build/'source-snapshot.zip') as archive:
        assert set(archive.namelist())==set(sources)
        for name,expected in sources.items():assert hashlib.sha256(archive.read(name)).hexdigest()==expected,name
    compiled={k:v for k,v in sources.items() if k.startswith(('src/','tests/','examples/')) or k in ('Cargo.toml','Cargo.lock','rust-toolchain.toml')}
    assert len(compiled)==204 and all(sha(root/name)==value for name,value in compiled.items())
    for directory in ('terrain-job','terrain-refined-job','terrain-ultrafine-job'):
        status=read(study/directory/'status.json');native=read(study/directory/'result.json')
        assert status['state']=='succeeded' and status['artifact_sha256']==sha(study/directory/'result.json')
        assert native['complete'] and native['schema']=='soil-fem-projected-mcc-stage-result/1'
        assert native['provenance']['solver_binary_sha256']==digest
    output.mkdir(exist_ok=False);(output/'solver').mkdir();shutil.copy2(binary,output/'solver/soil-fem-3d.exe')
    for name in ('source-snapshot.json','source-snapshot.zip'):shutil.copy2(build/name,output/name)
    (output/'logs').mkdir()
    for name in logs:shutil.copy2(runs/name,output/'logs'/name)
    (output/'receipts').mkdir()
    for name,path in [('lifecycle',lifecycle),('export',exported),('legacy-export',legacy_export),('legacy-solver',legacy_solver)]:
        shutil.copy2(path/'verification.json',output/'receipts'/(name+'.json'))
    shutil.copytree(assessment_dir,output/'refinement-assessment')
    verification=output/'verification-sources';verification.mkdir();freeze(verification)
    # Completed run folders are immutable evidence. Keep them by workspace
    # reference rather than duplicate large results and corrupted-input JSON.
    folders=[study,assessment_dir,lifecycle,exported,legacy_export,legacy_solver,
             runs/'projected-mcc-lifecycle-20261025-v1',runs/'projected-mcc-lifecycle-20261025-v2',
             runs/'projected-mcc-refinement-assessment-20261025-v1']
    files={p for folder in folders for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    files.update(p for p in output.rglob('*') if p.is_file())
    files.update(runs/name for name in logs)
    # Preserve failures and the actual frozen diagnostic/source versions.
    for name in ('projected-mcc-targeted-tests-20261025-v1.log','projected-mcc-targeted-tests-20261025-v2.log',
                 'projected-mcc-targeted-tests-20261025-v3.log','projected-mcc-targeted-tests-20261025-v4.log',
                 'projected-mcc-targeted-tests-20261025-v5.log','projected-mcc-cycle-reference-50-20261025-v1.log',
                 'projected-mcc-initial-reference-200-20261025-v1.log','projected-mcc-reference-400-driver-20261025-v1.log'):
        files.add(runs/name)
    report=dict(schema='soil-fem-projected-mcc-fixed-mesh-verification/1',passed=True,
        production_backend_accepted=False,general_mixed_terrain_accepted=False,
        goal_complete=False,p6_complete=False,spatial_convergence_verified=False,native_ui_connected=False,
        scope='Fixed actual CAD mesh, dry MCC P1-volume Newton; initial/load/unload/local excavation, frozen compatibility history, independently refined rate FEM, restart/rejection and dual-stress output',
        solver_sha256=digest,source_snapshot_sha256=sha(output/'source-snapshot.zip'),build_source_files=488,
        unchanged_compiled_source_files=204,verification_source_snapshot_sha256=sha(verification/'source-snapshot.zip'),
        release_rust_tests_passed=242,independent_reference_steps=assessment['independent_reference_steps'],
        fixed_mesh_integration_gates_passed=True,separate_process_midpoint_restarts=3,rejected_converged_trials=3,
        genuine_rejected_plastic_and_active_set_changes_verified=True,
        projected_export_invalid_input_refusals=14,legacy_export_invalid_input_refusals=15,
        legacy_mcc_mc_nonruntime_fields_exact=True,legacy_mcc_mc_export_files_byte_exact=True,
        finest_candidate=assessment['candidates'][-1],
        limitations=['General spatial convergence and comparison to an independently refined compatible formulation remain unverified',
                     'Water/member/JOINT/birth/material transitions and native viewer are not connected to this new formulation',
                     'Independent force-audited projected export currently supports explicit/side-roller kinematics',
                     'Some evidence is referenced in the workspace; this is not the P6 clean-environment distribution package',
                     'All other unfinished P0-P6 requirements and post-completion GitHub upload remain required'],
        sha256={p.relative_to(root).as_posix():sha(p) for p in sorted(files)})
    with (output/'final-validation.json').open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2,allow_nan=False)
    print(json.dumps({k:report[k] for k in ('passed','release_rust_tests_passed','independent_reference_steps',
        'fixed_mesh_integration_gates_passed','production_backend_accepted','p6_complete')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('assessment',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();record(a.study.resolve(),a.assessment.resolve(),a.output.resolve())
