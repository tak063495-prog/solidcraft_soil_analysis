"""Existing MCC accuracy gates on one CAD mesh, with both solvers refined.

Passing these gates cannot establish spatial convergence, coupled physics,
native integration or P6 completion. Coarse failures remain in the report.
"""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import numpy as np
import compare_projected_mcc as compare
import compare_mcc_rate_fem
import mixed_mcc_io
import fem_reference
from assess_projected_volume_diagnostic import fields


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def assess(study, output, coarse_steps=100, fine_steps=200):
    if not 2 <= coarse_steps < fine_steps:raise ValueError('ordered reference refinements required')
    def filename(n):return 'reference-cycle-100.json' if n==100 else f'reference-cycle-{n}-initial{n}.json'
    reference_coarse = study/filename(coarse_steps)
    reference_fine = study/filename(fine_steps)
    coarse, fine = read(reference_coarse), read(reference_fine)
    assert coarse['reference_steps'] == coarse_steps and fine['reference_steps'] == fine_steps
    assert coarse['initial_reference_sha256'] == sha(study/f'initial-{coarse_steps}.json')
    assert fine['initial_reference_sha256'] == sha(study/f'initial-{fine_steps}.json')
    for reference in (coarse, fine):
        for name, digest in reference['input_and_source_sha256'].items():assert sha(Path(name)) == digest, name
    output.mkdir(exist_ok=False)
    gates, refinements, candidates = [], [], []
    files = [reference_coarse, reference_fine, study/f'initial-{coarse_steps}.json', study/f'initial-{fine_steps}.json',
             Path(__file__), Path(compare.__file__), Path(compare_mcc_rate_fem.__file__),
             Path(mixed_mcc_io.__file__), Path(fem_reference.__file__)]

    def gate(name, value, limit, required=True):
        gates.append(dict(name=name, observed=float(value), limit=limit, required=required,
                          passed=bool(np.isfinite(value) and value < limit)))

    for a, b in zip(coarse['frames'], fine['frames']):
        assert a['name'] == b['name'] and a['active_elements'] == b['active_elements']
        difference = fields(a, b)
        refinements.append(dict(name=a['name'], fields=difference))
        for field, key, limit in (
            ('displacement_m', 'maximum_absolute', 1e-8), ('stress_kpa', 'relative_l2', 2.5e-6),
            ('pc_kpa', 'maximum_absolute', 2.5e-5), ('void_ratio', 'maximum_absolute', 1.25e-6),
            ('plastic_strain', 'maximum_absolute', 2.5e-7)):
            gate(f'independent{coarse_steps}_to_{fine_steps}:'+a['name']+':'+field, difference[field][key], limit)
        for label, frame in ((str(coarse_steps), a), (str(fine_steps), b)):
            gate('independent'+label+':'+a['name']+':force', frame['relative_work_conjugate_force_residual'], 1e-9)
    physics = None
    for label, directory, finest in (
        ('coarse', 'terrain-job', False), ('refined', 'terrain-refined-job', False),
        ('ultrafine', 'terrain-ultrafine-job', True)):
        job = study/directory;native_path = job/'result.json'
        native, manifest, status = map(read, (native_path, job/'job.json', job/'status.json'))
        assert status['state'] == 'succeeded' and status['artifact_sha256'] == sha(native_path)
        assert manifest['binding']['identity'] == native['provenance']
        for name, digest in manifest['binding']['files_sha256'].items():assert sha(job/name) == digest, name
        assert sha(job/'inputs/mesh.json') == sha(study/'inputs/mesh.json')
        assert sha(job/'inputs/ground.json') == sha(study/'inputs/ground.json')
        case = {k: v for k, v in native['case'].items() if k not in ('settings', 'verification_reject_once_stages')}
        case['kinematics']=native['case']['settings']['kinematics']
        case['stages'] = [{k: v for k, v in s.items() if k != 'initial_increment'} for s in case['stages']]
        if physics is None:physics = case
        else:assert case == physics, 'refinement changed physical case'
        target = output/('comparison-'+label+'.json')
        with contextlib.redirect_stdout(io.StringIO()):
            compare.compare(study/'inputs/mesh.json', reference_fine, native_path, target)
        comparison = read(target)
        for row, nf in zip(comparison['frames'], native['frames']):
            prefix = label+':'+row['name']+':'
            for name, limit in (('displacement_relative_all_nodes', 3e-5),
                                ('displacement_maximum_component_error_m', 1e-6),
                                ('independent_relative_force_residual', 1e-9)):
                gate(prefix+name, row[name], limit, finest or name == 'independent_relative_force_residual')
            for field, key, limit in (('stress', 'relative_all_points', 1e-5),
                ('pc', 'maximum_absolute_all_points', 1e-4), ('void_ratio', 'maximum_absolute_all_points', 5e-6),
                ('plastic_strain', 'maximum_absolute_all_points', 1e-6)):
                gate(prefix+field, row['fields'][field][key], limit, finest)
            tolerance = native['case']['settings']['relative_force_tolerance']*1.01
            gate(prefix+'native_force', nf['relative_force_residual'], tolerance)
            for name in ('relative_volume_equation_residual', 'relative_pressure_consistency_residual'):
                gate(prefix+name, nf['state']['mixed_mcc'][name], tolerance)
        final = native['frames'][-1]['monitors']['excavation_floor']
        rf = fine['frames'][-1]['monitors']['excavation_floor']
        candidates.append(dict(name=label, numerical_settings=native['case']['settings'],
            compared_frames=len(comparison['frames']), finest_candidate=finest,
            excavation_floor_increment_native_m=final['stage_increment_z_m'],
            excavation_floor_increment_independent_m=rf['mean_stage_increment_m'],
            comparison_file=target.name))
        files.extend([native_path, job/'job.json', job/'status.json', target])
    passed = all(row['passed'] for row in gates if row['required'])
    source_module = Path(__file__).with_name('assess_projected_volume_diagnostic.py')
    files.extend([source_module,Path(__file__).with_name('assess_mixed_mcc_resolution.py')])
    report = dict(schema='soil-fem-projected-mcc-integration-assessment/1',
        integration_gates_passed=passed, production_backend_accepted=False,
        spatial_convergence_verified=False, native_ui_connected=False, goal_complete=False, p6_complete=False,
        scope='One fixed actual CAD slope, dry MCC, P1-volume work-conjugate equilibrium; initialization/load/unload/local excavation',
        limits_origin='Unchanged MCC full-field and reference-refinement gates from assess_mixed_mcc_resolution.py; only the finest candidate is required for field accuracy and all coarse failures are retained',
        reference_refinement_includes_initialization=True, independent_reference_steps=[coarse_steps, fine_steps],
        candidates=candidates, independent_refinement=refinements, gates=gates,
        limitations='No mesh-convergence, water/member/JOINT coupling, oedometer export or native display acceptance',
        sha256={str(p.resolve()): sha(p) for p in files})
    with (output/'assessment.json').open('x', encoding='utf-8') as stream:json.dump(report, stream, indent=2)
    print(json.dumps(dict(integration_gates_passed=passed,
        failed_required=[g for g in gates if g['required'] and not g['passed']],
        retained_coarse_failures=sum(not g['passed'] and not g['required'] for g in gates),
        candidates=candidates)), flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('study', type=Path);p.add_argument('output', type=Path)
    p.add_argument('--coarse-steps',type=int,default=100);p.add_argument('--fine-steps',type=int,default=200)
    a=p.parse_args();assess(a.study.resolve(), a.output.resolve(), a.coarse_steps,a.fine_steps)
