"""Existing MCC accuracy/refinement gates for one explicitly bound CAD job.

A fine reference without a separately refined coarser reference can pass field
accuracy, but cannot earn integration or spatial-convergence acceptance.
"""
import argparse
import contextlib
import io
import json
from pathlib import Path
import numpy as np
import compare_projected_mcc as comparison
from assess_projected_volume_diagnostic import fields
import assess_projected_volume_diagnostic as field_module
import mcc_projected_cycle_diagnostic as cycle


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return comparison.sha(path)


def reference_bindings(path):
    report = read(path)
    assert report['schema'] == 'soil-fem-projected-mcc-cycle-diagnostic/1'
    bindings = report['input_and_source_sha256']
    for name, digest in bindings.items():
        assert sha(Path(name)) == digest, name
    initial = [Path(name) for name, digest in bindings.items() if digest == report['initial_reference_sha256']]
    assert len(initial) == 1
    initial_report = read(initial[0])
    assert initial_report['load_fraction'] == 1 and initial_report['reference_steps'] == report['reference_steps']
    for name, digest in initial_report['input_and_source_sha256'].items():
        assert sha(Path(name)) == digest, name
    return report, {Path(name).name: Path(name) for name in bindings}


def assess(job, reference, output, coarser_reference=None):
    status, manifest, native = [read(job/name) for name in ('status.json', 'job.json', 'result.json')]
    assert status['state'] == 'succeeded' and native['complete']
    assert status['artifact_sha256'] == sha(job/'result.json')
    assert manifest['binding']['identity'] == native['provenance']
    for name, digest in manifest['binding']['files_sha256'].items():
        assert sha(job/name) == digest, name
    fine, fine_paths = reference_bindings(reference)
    for name in ('mesh.json', 'ground.json'):
        assert sha(job/'inputs'/name) == sha(fine_paths[name]), name
    assert cycle.physical_case(read(job/'inputs/case.json')) == cycle.physical_case(read(fine_paths['case.json']))
    output.mkdir(parents=True, exist_ok=False)
    target = output/'comparison.json'
    with contextlib.redirect_stdout(io.StringIO()):
        comparison.compare(job/'inputs/mesh.json', reference, job/'result.json', target)
    compared = read(target)
    gates = []

    def gate(group, name, value, limit):
        gates.append(dict(group=group, name=name, observed=float(value), limit=limit,
                          passed=bool(np.isfinite(value) and value < limit)))

    assert len(compared['frames']) == len(native['frames']) == len(fine['frames'])
    for row, frame in zip(compared['frames'], native['frames']):
        name = row['name']
        for field, limit in (('displacement_relative_all_nodes', 3e-5),
                             ('displacement_maximum_component_error_m', 1e-6),
                             ('independent_relative_force_residual', 1e-9)):
            gate('field_accuracy', name+':'+field, row[field], limit)
        for field, key, limit in (('stress', 'relative_all_points', 1e-5), ('pc', 'maximum_absolute_all_points', 1e-4),
                                   ('void_ratio', 'maximum_absolute_all_points', 5e-6), ('plastic_strain', 'maximum_absolute_all_points', 1e-6)):
            gate('field_accuracy', name+':'+field, row['fields'][field][key], limit)
        tolerance = native['case']['settings']['relative_force_tolerance']*1.01
        gate('native_equations', name+':force', frame['relative_force_residual'], tolerance)
        for field in ('relative_volume_equation_residual', 'relative_pressure_consistency_residual'):
            gate('native_equations', name+':'+field, frame['state']['mixed_mcc'][field], tolerance)
    refinements = []
    if coarser_reference:
        coarse, paths = reference_bindings(coarser_reference)
        assert 2 <= coarse['reference_steps'] < fine['reference_steps']
        for name in ('mesh.json', 'ground.json', 'mesh-identity.json'):
            assert sha(paths[name]) == sha(fine_paths[name])
        assert cycle.physical_case(read(paths['case.json'])) == cycle.physical_case(read(fine_paths['case.json']))
        assert len(coarse['frames']) == len(fine['frames'])
        for a, b in zip(coarse['frames'], fine['frames']):
            assert a['name'] == b['name'] and a['active_elements'] == b['active_elements']
            difference = fields(a, b)
            refinements.append(dict(name=a['name'], fields=difference))
            for field, key, limit in (('displacement_m', 'maximum_absolute', 1e-8), ('stress_kpa', 'relative_l2', 2.5e-6),
                                      ('pc_kpa', 'maximum_absolute', 2.5e-5), ('void_ratio', 'maximum_absolute', 1.25e-6),
                                      ('plastic_strain', 'maximum_absolute', 2.5e-7)):
                gate('reference_refinement', a['name']+':'+field, difference[field][key], limit)
            for step, frame in ((coarse['reference_steps'], a), (fine['reference_steps'], b)):
                gate('reference_refinement', a['name']+f':force{step}', frame['relative_work_conjugate_force_residual'], 1e-9)
    accuracy = all(g['passed'] for g in gates if g['group'] in ('field_accuracy', 'native_equations'))
    refinement = all(g['passed'] for g in gates if g['group'] == 'reference_refinement') if coarser_reference else None
    evidence = [job/'status.json', job/'job.json', job/'result.json', reference, target,
                Path(__file__), Path(comparison.__file__), Path(field_module.__file__), Path(cycle.__file__)]
    if coarser_reference:
        evidence.append(coarser_reference)
    report = dict(schema='soil-fem-projected-mcc-single-job-assessment/1', field_accuracy_passed=accuracy,
                  reference_refinement_passed=refinement, integration_gates_passed=bool(accuracy and refinement is True),
                  reference_refinement_includes_initialization=bool(coarser_reference), reference_steps=fine['reference_steps'],
                  gates=gates, independent_refinement=refinements, production_backend_accepted=False,
                  spatial_convergence_verified=False, goal_complete=False, p6_complete=False,
                  limits_origin='Unchanged assess_projected_mcc_refinement.py and assess_mixed_mcc_resolution.py MCC full-field/refinement limits',
                  sha256={str(p): sha(p) for p in evidence},
                  limitations='One mesh and one native integration setting; further native refinement, spatial convergence and other P0-P6 requirements remain')
    with (output/'assessment.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(field_accuracy_passed=accuracy, reference_refinement_passed=refinement,
                          integration_gates_passed=report['integration_gates_passed'], failed=[g for g in gates if not g['passed']])), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('job', 'reference', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--coarser-reference', type=Path)
    args = parser.parse_args()
    assess(args.job.resolve(), args.reference.resolve(), args.output.resolve(),
           args.coarser_reference.resolve() if args.coarser_reference else None)
