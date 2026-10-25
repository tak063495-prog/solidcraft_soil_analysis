"""Seal the completed dry Tet10 plastic cycle without mixing build sources.

The numerical runs use the previously accepted production binary. Its primary
source archive is retained verbatim. Current verification sources, including
the newer independent midpoint implementation, are frozen separately and are
not presented as the build sources of that older production binary.
"""
import argparse
import json
from pathlib import Path
import shutil
import zipfile
import numpy as np
from compare_mcc_rate_fem import compare
from freeze_sources import freeze
from verify_acceptance_packet import sha,verify as verify_packet
from verify_mcc_cycle_lifecycle import no_elapsed


def read(path):return json.loads(path.read_text(encoding='utf-8'))


def record(out):
    root=Path(__file__).resolve().parents[1];run=root/'runs/mcc-rate-general-development-20261024-v1'
    parent=root/'runs/mcc-general-rate-acceptance-20261024-v1'
    lifecycle=root/'runs/mcc-cycle-lifecycle-20261025-v1';audit=root/'runs/mcc-cycle-audit-20261025-v1'
    parent_integrity=verify_packet(parent);parent_record=read(parent/'final-validation.json')
    solver=root/'target/release/soil-fem-3d.exe';binary_hash=sha(solver)
    assert binary_hash==parent_integrity['solver_sha256']==parent_record['solver_sha256']
    source=root/'runs/geostatic-test-20261024-003843-2c98f3fd/slope-pit/analysis/ground'
    mesh=read(source/'mesh.json');assert len(mesh['tetrahedra10'])==224 and len(mesh['nodes'])==459
    for name in ('mesh.json','ground.json'):assert sha(source/name)==sha(parent/'inputs'/name),name
    assessment=read(run/'cycle-global-refinement-assessment.json')
    assert assessment['passed'] and all(g['passed'] for g in assessment['gates'])
    for name,digest in assessment['sha256'].items():assert sha(Path(name))==digest,name
    references=[read(run/f'slope-cycle-reference-{n}.json') for n in (25,50,100)]
    candidate=read(run/'slope-cycle-rust-global-refined.json');assert candidate['complete']
    assert candidate['provenance']['solver_binary_sha256']==binary_hash
    comparison=compare(references[-1],candidate,mesh)
    assert comparison==assessment['comparisons'][-1]
    for i,row in enumerate(comparison['frames']):
        assert row['independent_plastic_points']==row['rust_plastic_points']==(7 if i==0 else 188)
        assert row['displacement_relative_all_nodes']<3e-5 and row['displacement_maximum_component_error_m']<1e-6
        assert row['fields']['plastic_strain']['maximum_absolute_all_points']<1e-6
        assert row['independent_relative_force_residual']<1e-9
    # Every independent input must match actual input bytes. Historical source
    # fingerprints may resolve into the production archive, never old inputs.
    with zipfile.ZipFile(parent/'source-snapshot.zip') as archive:
        production_sources=read(parent/'source-snapshot.json')
        for reference in references:
            for name,digest in reference['input_sha256'].items():
                path=Path(name);path=path if path.is_absolute() else root/path
                if path.is_file() and sha(path)==digest:continue
                relative=path.relative_to(root).as_posix()
                assert relative.startswith(('tools/','tests/')) and production_sources[relative]==digest,relative
                import hashlib
                assert hashlib.sha256(archive.read(relative)).hexdigest()==digest,relative
    restored=read(lifecycle/'rollback.json');control=read(lifecycle/'control.json');life=read(lifecycle/'verification.json')
    assert life['passed'] and len(life['separate_process_restarts'])==5
    for name,digest in life['input_and_source_sha256'].items():assert sha(Path(name))==digest,name
    for name,digest in life['evidence_sha256'].items():assert sha(lifecycle/name)==digest,name
    for key in ('frames','last_committed','reference_displacement_m'):assert restored[key]==control[key],key
    rejected=[a for a in restored['attempts'] if not a['accepted']]
    assert len(rejected)==3 and rejected[0]['rejected_trial_plastic_changed'] and rejected[-1]['rejected_trial_active_set_changed']
    assert all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    for count in (17,18,19,20,21):
        resumed=read(lifecycle/f'restart-{count}-result.json');assert resumed['complete']
        for key in ('frames','last_committed','reference_displacement_m'):assert resumed[key]==restored[key],(count,key)
        for key in ('geostatic','attempts'):assert no_elapsed(resumed[key])==no_elapsed(restored[key]),(count,key)
    audited=read(audit/'verification.json');assert audited['passed'] and audited['excavation_release_and_removed_history_checked']
    for name,digest in audited['input_and_source_sha256'].items():assert sha(Path(name))==digest,name
    for name,digest in audited['evidence_sha256'].items():assert sha(audit/name)==digest,name
    assert sha(audit/'result.json')==sha(run/'slope-cycle-rust-global-refined.json')
    removal=candidate['excavation_audits'][0]
    assert removal['removed_element_count']==21 and abs(removal['removed_volume_m3']-16)<1e-9
    assert all(f['independent_force_residual']<1e-8 and f['relative_moment_error']<1e-7 for f in audited['frames'])
    failed=read(run/'cycle-refinement-assessment.json')
    assert not failed['passed'] and any(not g['passed'] and g['limit']==1e-6 for g in failed['gates'])
    old_assessor=run/'cycle-refinement-assessment-v1-source.py'
    assert sha(old_assessor)==failed['sha256'][str(root/'tools/assess_mcc_cycle_refinement.py')]
    logs=['cached-integrated-tests.log','slope-cycle-reference-25.log','slope-cycle-reference-50.log','slope-cycle-reference-100.log',
          'slope-cycle-rust.log','slope-cycle-rust-finer.log','slope-cycle-rust-global-refined.log',
          'cycle-lifecycle.log','cycle-refinement-assessment.log','cycle-global-refinement-assessment.log']
    for name in logs:assert (run/name).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    assert 'Ran 9 tests' in (run/'cached-integrated-tests.log').read_text()
    assert (root/'runs/mcc-cycle-audit-20261025-v1.log').read_text().rstrip().endswith('Exit code: 0')
    maximum=max(f['maximum_absolute_engineering_strain_from_seed'] for f in candidate['frames'])
    assert maximum<.025
    out.mkdir(exist_ok=False);(out/'solver').mkdir();shutil.copyfile(solver,out/'solver/soil-fem-3d.exe')
    for name in ('source-snapshot.zip','source-snapshot.json'):shutil.copyfile(parent/name,out/name)
    shutil.copyfile(parent/'final-validation.json',out/'parent-final-validation.json')
    shutil.copytree(parent/'inputs',out/'inputs');(out/'data').mkdir()
    raw=[f'slope-cycle-reference-{n}.json' for n in (25,50,100)]
    raw+=['slope-cycle-rust.json','slope-cycle-rust-finer.json','slope-cycle-rust-global-refined.json',
          'cycle-refinement-assessment.json','cycle-global-refinement-assessment.json','cycle-refinement-assessment-v1-source.py']
    for name in raw+logs:shutil.copyfile(run/name,out/'data'/name)
    for name in ('slope-cycle-case.json','slope-cycle-finer-case.json','slope-cycle-global-refined-case.json','mesh-identity.json'):
        shutil.copyfile(run/name,out/'inputs'/name)
    shutil.copytree(lifecycle,out/'lifecycle');shutil.copytree(audit,out/'force-and-vtu-audit')
    shutil.copyfile(root/'runs/mcc-cycle-audit-20261025-v1.log',out/'force-and-vtu-audit.log')
    verification=out/'verification-sources';verification.mkdir();freeze(verification)
    final=comparison['frames'][-1];floor=final['monitors']['excavation_floor']
    report=dict(schema='soil-fem-mcc-tet10-cycle-final-validation/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Dry small-strain planar SolidCraft slope: plastic initial equilibrium, 5 kPa boundary load/unload and local excavation; all 224 Tet10/6048 GP',
        solver_sha256=binary_hash,source_snapshot_sha256=sha(out/'source-snapshot.zip'),
        production_source_origin='Verbatim previously accepted production archive; matches the retained production binary',
        parent_manifest_sha256=sha(parent/'final-validation.json'),
        verification_source_snapshot_sha256=sha(verification/'source-snapshot.zip'),
        verification_source_scope='Current verification workspace, including newer independent reference; not the production build source',
        reference_counts=[25,50,100],plastic_points_after_loading=188,
        reference_refinement=assessment['reference_refinement'],maximum_engineering_strain_from_seed=maximum,
        final_displacement_relative_error=final['displacement_relative_all_nodes'],
        final_displacement_maximum_component_error_m=final['displacement_maximum_component_error_m'],
        maximum_plastic_strain_component_error=max(f['fields']['plastic_strain']['maximum_absolute_all_points'] for f in comparison['frames']),
        floor_excavation_increment_independent_m=floor['independent_stage_increment_mean_m'],
        floor_excavation_increment_rust_m=floor['rust_stage_increment_mean_m'],
        converged_trial_rejections=3,separate_process_restarts=5,force_release_and_raw_vtu_audit_passed=True,
        limitations=['Numerical benchmark with about 2.10 percent peak engineering strain; no site calibration or large-deformation claim',
                     'No spatial/domain convergence, arbitrary curved plastic terrain, nonzero displacement paths or other material/coupling certification',
                     'Mixed elasticity and the current working-tree binary have separate evidence; all remaining P0-P6 requirements remain'],
        sha256={p.relative_to(root).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'final-validation.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:report[k] for k in ('passed','goal_complete','final_displacement_relative_error','maximum_plastic_strain_component_error','separate_process_restarts')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();record(a.output.resolve())
