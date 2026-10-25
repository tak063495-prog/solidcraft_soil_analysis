"""Seal the completed Tet10 slope equilibrium/excavation scope, preserving pending work."""
import argparse,json,shutil
from pathlib import Path
import numpy as np
import fem_reference as ref
from compare_mcc_rate_fem import compare
from freeze_sources import freeze
from record_mcc_rate_fem import no_elapsed, read
from verify_acceptance_packet import sha


def record(run,out):
    root=Path(__file__).resolve().parents[1]
    source=root/'runs/geostatic-test-20261024-003843-2c98f3fd/slope-pit/analysis/ground'
    mesh=read(source/'mesh.json');ground=read(source/'ground.json');solver=root/'target/release/soil-fem-3d.exe'
    parent=read(source.parents[2]/'final-validation.json')
    parent_hashes={n.replace('\\','/'):h for n,h in parent['sha256'].items()}
    names=['mesh.json','ground.json','bodies.json','body-0001.step','body-0002.step','mesh.provenance.json']
    for name in ('mesh.json','ground.json'):
        p=source/name;assert sha(p)==parent_hashes[p.relative_to(root).as_posix()],name
    accepted_path=source.parent/'accepted-stage-run.json'
    assert sha(accepted_path)==parent_hashes[accepted_path.relative_to(root).as_posix()]
    accepted={str(Path(n)):h for n,h in read(accepted_path)['inputs_and_results_sha256'].items()}
    for name in ('mesh.json','ground.json'):assert sha(source/name)==accepted[str(source/name)]
    bodies=read(source/'bodies.json');provenance=read(source/'mesh.provenance.json')
    assert provenance['cad_manifest_sha256']==sha(source/'bodies.json')
    for body in bodies['regions']+bodies['excavation_tools']:
        assert Path(body['step']).name==body['step']
        assert sha(source/body['step'])==body['source_sha256']
    by_body={b['body']:b for b in bodies['regions']}
    for region in mesh['regions']:assert region['source_sha256']==by_body[region['body']]['source_sha256']
    assert len(mesh['tetrahedra10'])==224 and len(mesh['nodes'])==459
    references=[read(run/f'slope-reference-{n}.json') for n in (25,50,100)]
    coarse=read(run/'slope-rust-fine.json');fine=read(run/'slope-rust-finer.json')
    comparison=compare(references[-1],fine,mesh);previous=compare(references[-1],coarse,mesh)
    assert comparison['integration_layout_checked_independently']
    source_hash=sha(root/'tools/mcc_rate_fem.py')
    for reference in references:assert reference['input_sha256'][str(root/'tools/mcc_rate_fem.py')]==source_hash
    for i,row in enumerate(comparison['frames']):
        assert row['independent_plastic_points']==row['rust_plastic_points']==7
        assert row['displacement_relative_all_nodes']<2e-5
        assert row['displacement_maximum_component_error_m']<3e-7
        assert row['display_displacement_relative_all_nodes']<2e-5
        assert row['independent_relative_force_residual']<1e-9
        assert row['displacement_relative_all_nodes']<previous['frames'][i]['displacement_relative_all_nodes']
        assert row['fields']['stress']['relative_all_points']<1e-5
        assert row['fields']['pc']['maximum_absolute_all_points']<1e-5
        assert row['fields']['void_ratio']['maximum_absolute_all_points']<5e-7
        assert row['fields']['plastic_strain']['maximum_absolute_all_points']<2e-7
    removed=mesh['excavation_sets']['Pit'];volume=sum(float(ref.volume(mesh,e)['weight'].sum()) for e in removed)
    assert len(removed)==21 and abs(volume-16)<1e-9
    a=np.asarray(references[-1]['frames'][0]['material_state']);b=np.asarray(references[-1]['frames'][1]['material_state'])
    assert np.array_equal(a[removed],b[removed])
    for e in removed:
        for q in range(27):assert fine['frames'][0]['state']['material_states'][e*27+q]==fine['frames'][1]['state']['material_states'][e*27+q]
    refinement=[]
    for i in range(2):
        arrays=[np.asarray(d['frames'][i]['displacement_m']) for d in references]
        changes=[float(np.linalg.norm(a-b)/np.linalg.norm(b)) for a,b in zip(arrays,arrays[1:])]
        assert changes[1]<2e-6 and changes[1]<.4*changes[0]
        refinement.append(dict(frame=comparison['frames'][i]['name'],relative_differences=changes))
    rollback=read(run/'slope-rollback.json');control=read(run/'slope-control.json')
    assert rollback['complete'] and control['complete']
    for key in ('frames','last_committed','reference_displacement_m'):assert rollback[key]==control[key]
    rejected=[a for a in rollback['geostatic']['attempts']+rollback['attempts'] if not a['accepted']]
    assert len(rejected)==2 and rejected[0]['rejected_trial_plastic_changed']
    assert all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    for n in (1,2,3):
        resumed=read(run/f'general-restart-{n}-result.json');assert resumed['complete']
        for key in ('frames','last_committed','reference_displacement_m'):assert resumed[key]==rollback[key]
        for key in ('geostatic','attempts'):assert no_elapsed(resumed[key])==no_elapsed(rollback[key])
    logs=['physics-tests-v2','identity','flat-elastic-25','slope-reference-25','slope-reference-50','slope-reference-100',
          'slope-rust-fine','slope-rust-finer','slope-rollback','slope-control','general-restart-verification']
    logs += [f'general-restart-{n}-{phase}' for n in (1,2,3) for phase in ('save','resume')]
    for name in logs:assert (run/(name+'.log')).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    assert 'Ran 8 tests' in (run/'physics-tests-v2.log').read_text()
    maximum=references[-1]['geostatic']['maximum_engineering_strain']
    assert maximum<ground['initialization']['geostatic']['maximum_initial_engineering_strain']
    out.mkdir(exist_ok=False);(out/'solver').mkdir();shutil.copyfile(solver,out/'solver/soil-fem-3d.exe')
    inputs=out/'inputs';inputs.mkdir();data=out/'data';data.mkdir()
    for name in names:shutil.copyfile(source/name,inputs/name)
    shutil.copyfile(accepted_path,inputs/'parent-accepted-stage-run.json')
    for name in ('slope-fine-case.json','slope-finer-case.json','slope-rollback-ground.json','slope-control-ground.json','slope-rollback-case.json','slope-control-case.json','mesh-identity.json'):
        shutil.copyfile(run/name,inputs/name)
    shutil.copyfile(root/'examples/stages-geostatic-slope-pit.json',inputs/'reference-case.json')
    raw=[f'slope-reference-{n}.json' for n in (25,50,100)]+['slope-rust-fine.json','slope-rust-finer.json','slope-rollback.json','slope-control.json']
    raw += [f'general-restart-{n}-{kind}.json' for n in (1,2,3) for kind in ('checkpoint','result')]
    raw += ['general-restart-verification.json','slope-rollback-verification.json','flat-elastic-25.json']
    for name in raw+[n+'.log' for n in logs]:shutil.copyfile(run/name,data/name)
    shutil.copyfile(root/'runs/check-mcc-general-restart-20261024.py',data/'check-restart.py')
    with (out/'comparison.json').open('x') as f:json.dump(comparison,f,indent=2)
    freeze(out)
    report=dict(schema='soil-fem-mcc-general-rate-final-validation/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Dry planar SolidCraft slope, all 224 Tet10/6048 GP, explicit point seeds, plastic initial equilibrium and local excavation; base case only',
        solver_sha256=sha(solver),source_snapshot_sha256=sha(out/'source-snapshot.zip'),reference_source_sha256=source_hash,
        maximum_initial_engineering_strain=maximum,removed_volume_m3=volume,plastic_points=7,reference_refinement=refinement,
        all_integration_positions_weights_and_owners_verified=True,rollback_and_three_process_restarts_exact=True,
        display_displacement_relative_error=comparison['frames'][-1]['display_displacement_relative_all_nodes'],
        pending=['5 kPa load/unload/excavation with additional plasticity is still running and is not included in this acceptance',
                 'arbitrary curved terrain, nonzero prescribed displacements, MC/IDH and coupled plastic displacement',
                 'spatial/domain convergence, post-localization objectivity, actual 2D comparison and all remaining P6 requirements'],
        sha256={p.relative_to(root).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'final-validation.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(dict(passed=True,goal_complete=False,all_points=6048,display_displacement_relative_error=report['display_displacement_relative_error'])))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();record(a.run.resolve(),a.output.resolve())
