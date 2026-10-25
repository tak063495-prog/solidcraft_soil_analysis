"""Bind fresh construction integration and the complete prior acceptance suite."""
import argparse
import json
from pathlib import Path
from freeze_sources import manifest
from verify_acceptance_packet import verify,sha

def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def record(run):
    root=Path(__file__).resolve().parents[1]
    assert read(run/'source-snapshot.json')==manifest(),'sources changed during acceptance'
    parent=verify(run/'baseline');old=read(run/'baseline/final-validation.json');cases=read(run/'cases/verification.json')
    digest=sha(run/'solver/soil-fem-3d.exe');assert digest==parent['solver_sha256']==cases['solver_sha256']==sha(root/'target/release/soil-fem-3d.exe')
    assert old['rust_tests']>=185 and old['coupled_separate_process_restarts']==29 and old['exact_legacy_cases']==6
    for name in ('baseline','projection-reference','cases'):assert (run/f'{name}.log').read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
    rust=(run/'baseline/baseline/core/rust-tests.log').read_text(encoding='utf-8')
    for name in ('saturated_flooded_cut_matches_exact_oedometer_unloading_and_mass_theory',
                 'domain_death_retains_remaining_heads_storage_fluxes_traces_and_archives_removed_water',
                 'rejected_domain_candidate_restores_history_and_clean_release_path_is_exact',
                 'construction_checkpoints_restore_nonzero_holding_domain_archive_and_next_physical_time',
                 'missing_cut_boundaries_overlap_reactivation_and_corrupt_origin_are_rejected',
                 'compressible_water_mass_removed_at_fixed_time_includes_positive_storage',
                 'two_local_removals_preserve_separate_archives_and_cumulative_water_budget',
                 'genuine_post_cut_hydraulic_exhaustion_preserves_the_last_accepted_domain_and_mass_archive',
                 'unchanged_domain_reuses_geometry_but_applies_each_new_period_boundary',
                 'curved_triangle_rain_projection_crossing_zero_normal_matches_exact_integral',
                 'planar_projection_handles_vertical_downward_and_scaled_shifted_geometry'):
        assert name+' ... ok' in rust,name
    assert cases['passed'] and cases['cad_exercised'] and len(cases['restarts'])==28
    assert all(row['all_states_frames_attempts_and_archives_exact'] for row in cases['restarts'])
    assert cases['cutback_clean_release_exact'] and cases['omitted_primary_exact'] and cases['invalid_future_boundary_rejected']
    assert cases['existing_output_preserved'] and cases['genuine_post_cut_exhaustion_preserves_committed_domain_and_archive']
    metrics={}
    for name in ('mc-tet4','mc-tet10','mcc-tet4','mcc-tet10','cutback','clean-release','omitted','cad-mc','cad-mcc'):
        folder=run/f'cases/{name}';data=read(folder/'result.json');assert data['complete'] and data['solver_binary_sha256']==digest
        assert read(folder/'verification.json')['passed'] and read(folder/'results/export-verification.json')['all_vtu_and_pvd_roundtrip']=='passed'
        if name.startswith('cad-'):
            assert read(folder/'cad-input-link.json')['mesh_contents_equal']
            mesh=read(folder/'mesh.json');state=data['last_committed'];initial=data['initialization']['last_committed']
            floor=mesh['node_sets']['Pit:floor'];surface=mesh['node_sets']['Pit:surrounding_surface']
            delta=[u[2]-v[2] for u,v in zip(state['mechanical']['displacement_m'],initial['mechanical']['displacement_m'])]
            metrics[name]=dict(nodes=len(mesh['nodes']),tet10=len(mesh['tetrahedra10']),removed_water_m3=data['excavations'][0]['removal']['water_removed_m3'],
                floor_nodal_mean_increment_mm=sum(delta[n] for n in floor)/len(floor)*1000,
                surrounding_surface_minimum_increment_mm=min(delta[n] for n in surface)*1000)
    files=[p for p in run.rglob('*') if p.is_file()]
    result=dict(schema='soil-fem-one-way-construction-final-validation/1',passed=True,solver_sha256=digest,source_snapshot_sha256=sha(run/'source-snapshot.zip'),
        rust_tests=old['rust_tests'],cad_tests=5,joint_split_tests=2,construction_separate_process_restarts=28,fixed_domain_coupled_restarts=29,standalone_hydraulic_restarts=8,exact_legacy_cases=6,
        parent_manifest_sha256=parent['manifest_sha256'],independent_elastic_comparisons=cases['elastic_comparisons'],cad_examples=metrics,
        maximum_independent_mechanical_residual=cases['maximum_independent_mechanical_residual'],maximum_global_water_budget_error_m3=cases['maximum_water_budget_error_m3'],
        scope='CAD DEATH with exact remaining water primary-state/trace transfer, archived removed water and GP history, explicit cut water/exterior policies, mechanical stress release at fixed time, physical-time paired flow/mechanics, failure rollback, separate-clock exports and exact restart. Independent elastic 3D PDE given water heads and saturated flooded-cut theory. No BIRTH/MATSET water policy, wet structures/JOINT, general independent plastic 3D, pond dynamics, or full P6 completion.',
        sha256={p.relative_to(root).as_posix():sha(p) for p in sorted(files)})
    with (run/'final-validation.json').open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k not in ('sha256','independent_elastic_comparisons')},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);record(p.parse_args().run.resolve())
