"""Seal two independently refined dry MCC/Tet4 3D benchmarks.

This is a limited acceptance, never completion of P6 or all plastic models.
Acceptance checks recompute comparisons from raw states, not green flags.
"""
import argparse,hashlib,json,shutil,sys
from pathlib import Path
import numpy as np
from compare_mcc_rate_fem import compare
import fem_reference as ref
from freeze_sources import freeze
from verify_acceptance_packet import sha


def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def relative(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))
def no_elapsed(value):
    if isinstance(value,dict):return {k:no_elapsed(v) for k,v in value.items() if k!='elapsed_seconds'}
    if isinstance(value,list):return [no_elapsed(v) for v in value]
    return value


def record(run,out):
    root=Path(__file__).resolve().parents[1];mesh_path=root/'runs/local-verified-20261023/ground/mesh.json'
    ground_path=mesh_path.parent/'ground.json';mesh=read(mesh_path);solver=root/'target/release/soil-fem-3d.exe'
    assert not mesh.get('tetrahedra10') and len(mesh['tetrahedra'])==784
    cad_acceptance=read(mesh_path.parent/'accepted-run.json')
    for name,digest in cad_acceptance['files_sha256'].items():assert sha(mesh_path.parent/name)==digest,name
    provenance=read(mesh_path.parent/'mesh.provenance.json')
    assert provenance['cad_manifest_sha256']==sha(mesh_path.parent/'bodies.json')
    assert provenance['ground_sha256']==sha(ground_path)
    removed_volume=sum(float(ref.volume(mesh,int(e))['weight'].sum()) for e in mesh['excavation_sets']['Pit'])
    assert abs(removed_volume-32.)<1e-9
    references=[read(run/f'plastic-{n}.json') for n in (25,50,100)]
    natives=[read(run/f'rust-{n}.json') for n in ('fine','finer','finest','global-refined')]
    # The actual reference code remains bound to all completed executions.
    source_hash=sha(root/'tools/mcc_rate_fem.py')
    for reference in references+[read(run/'near-nc-reference-25.json'),read(run/'near-nc-reference-100.json')]:
        assert reference['input_sha256'][str(root/'tools/mcc_rate_fem.py')]==source_hash
    main=compare(references[-1],natives[-1],mesh)
    near=compare(read(run/'near-nc-reference-100.json'),read(run/'near-nc-rust-finer.json'),mesh)
    for comparison in (main,near):
        for row in comparison['frames'][1:]:
            assert row['displacement_relative_all_nodes']<1e-3,(row['name'],'relative displacement')
            assert row['displacement_maximum_component_error_m']<1e-6,(row['name'],'absolute displacement')
            assert row['independent_relative_force_residual']<1e-9
            for field in ('stress','pc'):assert row['fields'][field]['relative_all_points']<1e-5,(row['name'],field)
            assert row['fields']['void_ratio']['maximum_absolute_all_points']<1e-6
            assert row['fields']['plastic_strain']['maximum_absolute_all_points']<5e-7
        floor=comparison['frames'][-1]['monitors']['excavation_floor']
        surface=comparison['frames'][-1]['monitors']['surrounding_surface']
        assert floor['independent_stage_increment_mean_m']>0 and floor['rust_stage_increment_mean_m']>0
        assert surface['independent_minimum_from_seed_m']<0 and surface['rust_minimum_from_seed_m']<0
    assert main['frames'][1]['independent_plastic_points']==main['frames'][1]['rust_plastic_points']==132
    assert near['frames'][-1]['independent_plastic_points']==near['frames'][-1]['rust_plastic_points']==42
    refinement=[]
    for i in range(1,4):
        u=[np.asarray(r['frames'][i]['displacement_m']) for r in references]
        first=relative(u[0],u[1]);second=relative(u[1],u[2])
        assert second<4e-5 and second<.35*first,(i,first,second)
        errors=[relative(np.asarray(n['frames'][i]['state']['displacement_m']),u[-1]) for n in natives]
        assert all(b<a for a,b in zip(errors,errors[1:])),(i,errors)
        refinement.append(dict(stage=main['frames'][i]['name'],reference_counts=[25,50,100],reference_relative_differences=[first,second],
                               rust_relative_displacement_errors=errors))
    near_u=[np.asarray(read(run/f'near-nc-reference-{n}.json')['frames'][-1]['displacement_m']) for n in (25,100)]
    near_refinement=relative(*near_u);assert near_refinement<1e-4
    near_errors=[relative(np.asarray(read(run/f'near-nc-rust{suffix}.json')['frames'][-1]['state']['displacement_m']),near_u[-1]) for suffix in ('','-finer')]
    assert near_errors[1]<near_errors[0]
    maximum_strain=0.
    for native in (natives[-1],read(run/'near-nc-rust-finer.json')):
        assert native['complete']
        for frame in native['frames']:
            state=frame['state'];u=np.asarray(state['displacement_m'])
            for e in np.flatnonzero(state['active_elements']):
                q=ref.volume(mesh,int(e));maximum_strain=max(maximum_strain,float(np.max(np.abs(q['b']@u[q['nodes']].reshape(-1)))))
    assert maximum_strain<.01
    rejected=read(run/'rollback.json');control=read(run/'rollback-control.json')
    for key in ('last_committed','frames'):assert rejected[key]==control[key]
    trials=[a for a in rejected['attempts'] if not a['accepted']]
    assert len(trials)==2 and all(a['committed_state_unchanged_on_rejection'] for a in trials)
    assert trials[0]['rejected_trial_plastic_changed'] and trials[1]['rejected_trial_material_changed']
    for count in (2,6,11):
        resumed=read(run/f'restart-{count}-result.json');assert resumed['complete']
        for key in ('last_committed','frames'):assert resumed[key]==rejected[key]
        assert no_elapsed(resumed['attempts'])==no_elapsed(rejected['attempts'])
    log_names=['material-rate-tests-v2','plastic-25','plastic-50','plastic-100','rust-fine','rust-finer','rust-finest','rust-global-refined',
               'near-nc-reference-25','near-nc-reference-100','near-nc-rust','near-nc-rust-finer','rollback','rollback-control','restart-verification-v2']
    log_names += [f'restart-{n}-{kind}' for n in (2,6,11) for kind in ('save','resume')]
    for name in log_names:assert (run/(name+'.log')).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    assert 'Ran 4 tests' in (run/'material-rate-tests-v2.log').read_text()
    out.mkdir(exist_ok=False);(out/'solver').mkdir();shutil.copyfile(solver,out/'solver/soil-fem-3d.exe')
    data=out/'data';data.mkdir();inputs=out/'inputs';inputs.mkdir()
    shutil.copyfile(mesh_path,inputs/'mesh.json');shutil.copyfile(ground_path,inputs/'ground.json')
    cad=inputs/'cad';cad.mkdir()
    for name in cad_acceptance['files_sha256']:shutil.copyfile(mesh_path.parent/name,cad/name)
    shutil.copyfile(mesh_path.parent/'accepted-run.json',cad/'accepted-run.json')
    for path in run.glob('*case.json'):shutil.copyfile(path,inputs/path.name)
    shutil.copyfile(run/'near-nc-ground.json',inputs/'near-nc-ground.json')
    raw_names=[f'plastic-{n}.json' for n in (25,50,100)]+[f'rust-{n}.json' for n in ('fine','finer','finest','global-refined')]
    raw_names += ['near-nc-reference-25.json','near-nc-reference-100.json','near-nc-rust.json','near-nc-rust-finer.json','rollback.json','rollback-control.json']
    raw_names += [f'restart-{n}-{kind}.json' for n in (2,6,11) for kind in ('checkpoint','result')]
    for name in raw_names+[n+'.log' for n in log_names]:shutil.copyfile(run/name,data/name)
    shutil.copyfile(root/'runs/check-mcc-rate-restart-20261024.py',data/'check-restart.py')
    for name,value in (('main-comparison.json',main),('near-nc-comparison.json',near)):
        with (out/name).open('x') as f:json.dump(value,f,indent=2)
    freeze(out)
    report=dict(schema='soil-fem-mcc-rate-fem-final-validation/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Two dry MCC/Tet4 full 3D cases on the same accepted SolidCraft mesh; plastic preload/unload/local excavation and new plasticity during excavation',
        solver_sha256=sha(solver),source_snapshot_sha256=sha(out/'source-snapshot.zip'),reference_source_sha256=source_hash,
        maximum_engineering_strain=maximum_strain,reference_and_rust_refinement=refinement,
        independently_integrated_removed_volume_m3=removed_volume,
        near_nc_reference_relative_refinement=near_refinement,near_nc_rust_relative_errors=near_errors,
        plastic_points_preload=132,plastic_points_excavation=42,rollback_all_frames_and_history_exact=True,
        separate_process_restart_points=[2,6,11],restart_timing_exclusion='elapsed_seconds only; all other attempt values exact',
        limits=['Tet10, general terrain and boundary conditions, MC/IDH and coupled plastic displacement still require independent acceptance',
          'No spatial/domain convergence claim: identical fixed Tet4 mesh compared',
          'No shear-band objectivity, post-localization uniqueness, wet/structural/SRM validation or full P6 completion claim'],
        sha256={p.relative_to(root).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'final-validation.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(dict(passed=True,goal_complete=False,maximum_engineering_strain=maximum_strain,
                         finest_stage_errors=[f['displacement_relative_all_nodes'] for f in main['frames']],
                         near_nc_error=near['frames'][-1]['displacement_relative_all_nodes'])))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();record(a.run.resolve(),a.output.resolve())
