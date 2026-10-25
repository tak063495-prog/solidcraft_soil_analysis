"""Full-field comparison of the explicitly identified projected MCC formulation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import compare_mcc_rate_fem as comparison
import fem_reference
import mixed_mcc_io


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def metric(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        raise ValueError('projected field dimensions differ')
    return dict(relative_l2=float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-30)),
                maximum_absolute=float(abs(a-b).max()))


def compare(mesh_file, reference_file, native_file, output):
    mesh, reference, native = map(read, (mesh_file, reference_file, native_file))
    if reference['schema'] != 'soil-fem-projected-mcc-cycle-diagnostic/1' or native['schema'] != mixed_mcc_io.PROJECTED_SCHEMA:
        raise ValueError('explicit projected MCC reference and completed native result required')
    mixed_mcc_io.validate(mesh, native)
    for name, digest in reference['input_and_source_sha256'].items():
        assert sha(Path(name)) == digest, name
    # The shared comparator consumes force residuals but does not interpret
    # their definition. Give it the explicitly labelled work-conjugate value.
    view = {**reference, 'frames': [{**f, 'relative_free_force_residual': f['relative_work_conjugate_force_residual']}
                                   for f in reference['frames']]}
    report = comparison.compare(view, native, mesh)
    for row, rf, nf in zip(report['frames'], reference['frames'], native['frames']):
        mixed = nf['state']['mixed_mcc'];projected = mixed['projected_volume']
        for label, actual, expected, field in (
            ('primary_mean_pressure_kpa', mixed['pressure_dofs'], rf['pressure_dofs'], 'compression_positive_kpa'),
            ('primary_volume_strain', projected['volume_dofs'], rf['volume_dofs'], 'compression_positive'),
            ('frozen_compatibility_offset_m3', projected['compatibility_offsets'], rf['compatibility_offsets'], 'integral_m3')):
            a = {(d['node_index'], d['material_id']): d[field] for d in actual}
            b = {(d['node_index'], d['material_id']): d[field] for d in expected}
            assert len(a) == len(actual) and len(b) == len(expected) and a.keys() == b.keys()
            row[label] = metric([a[k] for k in sorted(a)], [b[k] for k in sorted(b)])
    assert len(native['excavation_audits']) == len(reference['release_audits'])
    releases = []
    for actual, expected in zip(native['excavation_audits'], reference['release_audits']):
        assert actual['stage'] == expected['stage']
        releases.append(dict(stage=actual['stage'], holding_force_kn=metric(actual['initial_holding_load_kn'],
            np.asarray(expected['holding_nodal_load_kn']).reshape(-1)),
            native_zero_release_relative_force_residual=actual['zero_release_relative_force_residual']))
    files = [mesh_file, reference_file, native_file, Path(__file__), Path(comparison.__file__),
             Path(fem_reference.__file__), Path(mixed_mcc_io.__file__)]
    report.update(schema='soil-fem-projected-mcc-independent-comparison/1',
        development_comparison_only=True, production_backend_accepted=False, goal_complete=False, p6_complete=False,
        formulation=mixed_mcc_io.PROJECTED_FORMULATION,
        equilibrium='Work-conjugate material deviator + projected mean effective stress; raw material stress is retained separately',
        reference_steps=reference['reference_steps'], native_provenance=native['provenance'],
        release_audits=releases, independently_validated_primary_fields_and_force=True,
        limitations='Matching formulation on one mesh; both integration refinement and spatial convergence remain required',
        sha256={str(p.resolve()): sha(p) for p in files})
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(frames=[dict(name=f['name'], displacement_relative=f['displacement_relative_all_nodes'],
        displacement_maximum_error_m=f['displacement_maximum_component_error_m'],
        pc_maximum_error_kpa=f['fields']['pc']['maximum_absolute_all_points'],
        plastic_strain_maximum_error=f['fields']['plastic_strain']['maximum_absolute_all_points']) for f in report['frames']],
        production_backend_accepted=False)), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ('mesh', 'reference', 'native', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    compare(*(getattr(a, k).resolve() for k in ('mesh', 'reference', 'native', 'output')))
