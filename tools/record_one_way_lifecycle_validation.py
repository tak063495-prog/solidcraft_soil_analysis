"""Seal the fresh lifecycle evidence and unchanged parent acceptance sources."""
import argparse
import json
from pathlib import Path
from freeze_sources import manifest
from verify_acceptance_packet import verify,sha

def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def record(run):
    root=Path(__file__).resolve().parents[1];assert read(run/'source-snapshot.json')==manifest(),'sources changed during acceptance'
    parent=verify(run/'baseline');old=read(run/'baseline/final-validation.json');cases=read(run/'cases/verification.json');digest=sha(run/'solver/soil-fem-3d.exe')
    assert cases['passed'] and cases['cad_exercised'] and old['rust_tests']>=197
    assert digest==parent['solver_sha256']==cases['solver_sha256']==sha(root/'target/release/soil-fem-3d.exe')
    assert len(cases['restarts'])==33 and all(c['all_states_frames_attempts_and_archives_exact'] for c in cases['restarts'])
    assert old['construction_separate_process_restarts']==28 and old['fixed_domain_coupled_restarts']==29 and old['standalone_hydraulic_restarts']==8 and old['exact_legacy_cases']==6
    for name in ('baseline','cases'):assert (run/f'{name}.log').read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
    for name in ('plastic_rollback_clean_path_exact','mcc_archive_survives_inactive_replacement_and_birth','invalid_future_boundary_rejected','genuine_infeasible_water_transfer_preserves_last_committed','omitted_primary_exact'):assert cases[name]
    for folder in (run/'cases').iterdir():
        if folder.is_dir() and (folder/'result.json').exists():
            assert read(folder/'result.json')['solver_binary_sha256']==digest
            assert read(folder/'verification.json')['passed'] and read(folder/'results/export-verification.json')['all_vtu_and_pvd_roundtrip']=='passed'
    assert read(run/'cases/cad-wrapper/run/results/export-verification.json')['all_vtu_and_pvd_roundtrip']=='passed'
    files=sorted(p for p in run.rglob('*') if p.is_file())
    result=dict(schema='soil-fem-one-way-lifecycle-final-validation/1',passed=True,solver_sha256=digest,source_snapshot_sha256=sha(run/'source-snapshot.zip'),parent_manifest_sha256=parent['manifest_sha256'],
        rust_tests=old['rust_tests'],cad_tests=5,joint_split_tests=2,lifecycle_separate_process_restarts=33,construction_separate_process_restarts=28,fixed_domain_coupled_restarts=29,standalone_hydraulic_restarts=8,exact_legacy_cases=6,
        profiles=len(cases['profiles']),saturated_birth_theory=cases['saturated_birth_theory'],independent_elastic_comparisons=cases['elastic_comparisons'],
        maximum_independent_mechanical_residual=max(v['maximum_relative_mechanical_residual'] for v in cases['profiles'].values()),maximum_independent_water_budget_error_m3=max(v['maximum_global_water_budget_error_m3'] for v in cases['profiles'].values()),
        scope='Initially inactive soil, genuine time-zero excavation, fresh MC BIRTH and explicit water head, active/inactive MATSET with preserve-head or preserve-water law transfer, archived old histories, paired rollback, separate-process restart, actual SolidCraft prepared mesh and verified exports. Independent incremental elastic FE given hydraulic heads and saturated birth theory. No wet JOINT/structures, stress-free MCC birth, full plastic 3D reference, pond dynamics or full P6 completion.',sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    with (run/'final-validation.json').open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k not in ('sha256','independent_elastic_comparisons')},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',required=True,type=Path);record(p.parse_args().run.resolve())
