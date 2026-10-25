"""Roundtrip both stress definitions, nodal volume and frozen domain offsets."""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import meshio
import numpy as np
import export_stages
import high_order_io
import mixed_mcc_io
import fem_reference


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(study, output):
    root = Path(__file__).resolve().parents[1]
    cases = [('terrain', study/'inputs/mesh.json', study/'terrain-job/result.json'),
             ('two-material', root/'examples/hardening-column-tet10-mesh.json', study/'two-material-result.json')]
    sources = [Path(__file__), Path(export_stages.__file__), Path(high_order_io.__file__), Path(mixed_mcc_io.__file__), Path(fem_reference.__file__)]
    inputs = [p for _, m, r in cases for p in (m, r)]+sources
    fingerprints = {str(p): sha(p) for p in inputs}
    output.mkdir(exist_ok=False)
    reports = []
    for label, mesh_file, result_file in cases:
        mesh, result = read(mesh_file), read(result_file)
        target = output/label
        export_stages.export(mesh_file, result_file, target)
        metadata = read(target/'export-verification.json')
        assert metadata['schema'] == 'soil-fem-stage-export/4'
        assert metadata['equilibrium_stress_field'] == 'work_conjugate_effective_stress_compression_positive_kpa'
        gaps, jumps, offsets = [], [], []
        for index, frame in enumerate(result['frames']):
            fields = mixed_mcc_io.frame_fields(mesh, result, frame)
            primary = meshio.read(target/f'mean-pressure-dofs-{index:03}.vtu')
            for name, expected in (
                ('mean_effective_stress_compression_positive_kpa', fields['primary']),
                ('cumulative_material_volume_strain_compression_positive', fields['projected_volume_dofs']),
                ('frozen_compatibility_offset_m3', fields['projected_compatibility_offsets'])):
                assert np.array_equal(primary.point_data[name], expected), name
            raw = meshio.read(target/f'integration-points-{index:03}.vtu')
            ids = raw.point_data['integration_point_id']
            work = fields['raw']['work_conjugate_effective_stress_compression_positive_kpa'][ids]
            material = np.asarray([s['effective_stress_compression_positive_kpa'] for s in frame['state']['material_states']])[ids]
            assert np.array_equal(raw.point_data['work_conjugate_effective_stress_compression_positive_kpa'], work)
            assert np.array_equal(raw.point_data['effective_stress_compression_positive_kpa'], material)
            gaps.append(float(abs(work-material).max()))
            offsets.append(float(abs(fields['projected_compatibility_offsets']).max()))
            by_node = {}
            for key, p in zip(fields['keys'], fields['primary']):by_node.setdefault(key[0], []).append(p)
            jumps.extend(max(v)-min(v) for v in by_node.values() if len(v)>1)
            with (target/f'material-history-{index:03}.csv').open(encoding='utf-8-sig', newline='') as stream:
                history = list(csv.DictReader(stream))
            components = ('xx', 'yy', 'zz', 'xy', 'yz', 'zx')
            for point, row in enumerate(history):
                values = [row['work_conjugate_effective_stress_'+c+'_kpa'] for c in components]
                if row['active'] == '0':
                    assert values == ['']*6
                else:
                    assert np.array_equal([float(v) for v in values], fields['raw']['work_conjugate_effective_stress_compression_positive_kpa'][point])
        if label == 'terrain':assert max(gaps)>1. and max(offsets)>1e-5
        else:assert max(jumps)>1.
        reports.append(dict(name=label, frames=len(result['frames']), both_stress_definitions_roundtrip=True,
            nodal_volume_and_compatibility_offset_roundtrip=True, maximum_stress_definition_difference_kpa=max(gaps),
            maximum_material_pressure_jump_kpa=max(jumps, default=0.), maximum_frozen_offset_m3=max(offsets)))
    original = read(cases[0][2]);refusals = []

    def refuse(name, change):
        value = copy.deepcopy(original);change(value)
        path = output/(name+'.json');target = output/(name+'-must-not-exist')
        with path.open('x', encoding='utf-8') as stream:json.dump(value, stream)
        try:export_stages.export(cases[0][1], path, target)
        except ValueError as error:
            assert not target.exists(), name
            refusals.append(dict(name=name, reason=str(error), no_output_created=True))
        else:raise AssertionError(name+' accepted')

    def projection(d):return d['frames'][1]['state']['mixed_mcc']['projected_volume']
    refuse('missing-projection', lambda d: d['frames'][1]['state']['mixed_mcc'].pop('projected_volume'))
    refuse('unknown-formulation', lambda d: projection(d).update(formulation='unknown'))
    refuse('missing-volume-node', lambda d: projection(d)['volume_dofs'].pop())
    refuse('duplicate-volume-node', lambda d: projection(d)['volume_dofs'].append(copy.deepcopy(projection(d)['volume_dofs'][0])))
    refuse('changed-volume-node', lambda d: projection(d)['volume_dofs'][0].update(compression_positive=.5))
    refuse('missing-offset', lambda d: projection(d)['compatibility_offsets'].pop())
    refuse('changed-offset', lambda d: projection(d)['compatibility_offsets'][0].update(integral_m3=.01))
    refuse('false-pressure-gap', lambda d: projection(d).update(maximum_gp_pressure_projection_difference_kpa=100.))
    refuse('false-force-residual', lambda d: d['frames'][1].update(relative_force_residual=.1))
    refuse('false-reaction', lambda d: d['frames'][1]['support_reaction_kn'].__setitem__(2, 0.))
    refuse('false-applied-force', lambda d: d['frames'][1]['applied_force_kn'].__setitem__(2, 0.))
    refuse('changed-pressure-node', lambda d: d['frames'][1]['state']['mixed_mcc']['pressure_dofs'][0].update(compression_positive_kpa=200.))
    def relabel(d):
        d.update(schema=mixed_mcc_io.SCHEMA);d['case'].update(schema='soil-fem-mixed-mcc-stages/1')
    refuse('relabelled-as-gp-volume', relabel)
    def removed(d):
        state=d['frames'][-1]['state'];i=next(i for i,p in enumerate(d['integration_points']) if not state['active_elements'][p['element_index']])
        state['material_states'][i]['plastic_strain_compression_positive'][0]+=.001
    refuse('changed-frozen-material-history', removed)
    assert fingerprints == {str(p): sha(p) for p in inputs}
    report = dict(schema='soil-fem-projected-mcc-export-verification/1', passed=True,
        production_backend_accepted=False, goal_complete=False, p6_complete=False,
        scope='Output semantics, independent force/reaction and volume audits, roundtrips and corrupted-input refusal; no spatial physics acceptance',
        cases=reports, rejected_inputs=refusals, input_and_source_sha256=fingerprints,
        evidence_sha256={str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file()})
    with (output/'verification.json').open('x', encoding='utf-8') as stream:json.dump(report, stream, indent=2)
    print(json.dumps(dict(passed=True, cases=reports, rejected_inputs=len(refusals))), flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('study', type=Path);p.add_argument('output', type=Path)
    a=p.parse_args();verify(a.study.resolve(), a.output.resolve())
