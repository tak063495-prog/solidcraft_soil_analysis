"""Dedicated schema-15 raw data writer. No mechanical acceptance is inferred.

Common material-field averaging is reused from the unchanged general writer;
existing exporters and historical verification dependencies remain untouched.
"""
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import meshio
import numpy as np
from high_order_io import average, groups, state_fields

def export(mesh_file, result_file, output):
    fingerprints = {str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (mesh_file, result_file)}
    mesh = json.loads(Path(mesh_file).read_text(encoding='utf-8-sig'))
    result = json.loads(Path(result_file).read_text(encoding='utf-8-sig'))
    if result.get('schema') != 'soil-fem-stage-result/15' or not result.get('complete'):
        raise ValueError('completed contact-BIRTH result/15 required')
    if any((f['state'].get('mixed_mcc') is not None for f in result['frames'])):
        raise ValueError('mixed primary state requires its own result schema')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    quadratic = bool(mesh.get('tetrahedra10'))
    kind = 'tetra10' if quadratic else 'tetra'
    xyz = np.array(mesh['nodes'])
    cells = np.array(mesh['tetrahedra10'] if quadratic else mesh['tetrahedra'])
    layout = result['integration_points']
    ids = groups(layout, len(cells))
    initial_reference = np.asarray(result.get('reference_displacement_m', np.zeros_like(xyz)))
    roots = [ET.Element('VTKFile', type='Collection', version='0.1', byte_order='LittleEndian') for _ in range(2)]
    collections = [ET.SubElement(r, 'Collection') for r in roots]
    rows = []
    for frame_index, frame in enumerate(result['frames']):
        state = frame['state']
        reference = np.asarray(state.get('display_reference_displacement_m') or initial_reference)
        active = np.asarray(state['active_elements'])
        used = np.unique(cells[active])
        mapping = np.full(len(xyz), -1)
        mapping[used] = np.arange(len(used))
        u = np.asarray(state['displacement_m'])
        old = np.asarray(result['frames'][max(0, frame_index - 1)]['state']['displacement_m'])
        if frame_index and result['case']['stages'][frame_index - 1]['kind'] == 'birth':
            old = np.asarray(next((a for a in result['transition_audits'] if a['stage'] == frame['name']))['prepared_state']['displacement_m'])
        raw = state_fields(state['material_states'])
        h = state.get('hydrostatic_water_level_m', result['ground']['initialization']['water_level_m'])
        gamma = result['ground']['initialization']['water_unit_weight_kn_m3']
        pore = np.zeros(len(layout)) if h is None else np.maximum(h - np.array([p['position_m'][2] for p in layout]), 0) * gamma
        if state.get('pore_pressure_kpa') is not None:
            pore = np.asarray(state['pore_pressure_kpa'])
        structures = set((result['case'].get('joint_system') or {}).get('structure_bodies', []))
        bodies = {r['id']: r['body'] for r in mesh['regions']}
        total = np.asarray([bodies[mesh['region_ids'][p['element_index']]] in structures for p in layout])
        if state.get('pore_pressure_kpa') is not None:
            assert np.all(pore[total] == 0.0)
        pore = pore.copy()
        pore[total] = 0.0
        raw['material_stress_reference'] = total.astype(np.uint8)
        raw['pore_pressure_kpa'] = pore
        raw['total_stress_compression_positive_kpa'] = raw['effective_stress_compression_positive_kpa'].copy()
        raw['total_stress_compression_positive_kpa'][:, :3] += pore[:, None]
        assert all((len(v) == len(layout) for v in raw.values()))
        cell = {'original_element_id': [np.flatnonzero(active)], 'material_id': [np.asarray(state.get('element_material_ids') or mesh['material_ids'])[active]], 'region_id': [np.asarray(mesh['region_ids'])[active]], 'integration_point_count': [np.array([len(g) for g in ids])[active]]}
        cell.update({'volume_average_' + name: [average(value, layout, ids)[active]] for name, value in raw.items()})
        filename = f'stage-{frame_index:03}.vtu'
        node_fields = {'original_node_id': used, 'displacement_m': (u - reference)[used], 'stage_displacement_increment_m': (u - old)[used]}
        node_fields.update(total_displacement_from_seed_m=u[used], initialization_displacement_m=initial_reference[used])
        node_fields['display_reference_displacement_m'] = reference[used]
        meshio.write_points_cells(output / filename, xyz[used], [(kind, mapping[cells[active]])], point_data=node_fields, cell_data=cell)
        qp = np.array([i for i, p in enumerate(layout) if active[p['element_index']]])
        positions = np.array([layout[i]['position_m'] for i in qp])
        point_u = np.array([np.array(layout[i]['shape_functions']) @ (u - reference)[cells[layout[i]['element_index']]] for i in qp])
        point_file = f'integration-points-{frame_index:03}.vtu'
        point_fields = {**{k: v[qp] for k, v in raw.items()}, 'integration_point_id': qp, 'original_element_id': np.array([layout[i]['element_index'] for i in qp]), 'integration_weight_m3': np.array([layout[i]['weight_m3'] for i in qp]), 'displacement_m': point_u}
        point_fields['material_id'] = np.array([(state.get('element_material_ids') or mesh['material_ids'])[layout[i]['element_index']] for i in qp])
        point_fields['total_displacement_from_seed_m'] = np.array([np.array(layout[i]['shape_functions']) @ u[cells[layout[i]['element_index']]] for i in qp])
        point_fields['initialization_displacement_m'] = np.array([np.array(layout[i]['shape_functions']) @ initial_reference[cells[layout[i]['element_index']]] for i in qp])
        meshio.write_points_cells(output / point_file, positions, [('vertex', np.arange(len(qp))[:, None])], point_data=point_fields)
        checked = meshio.read(output / filename)
        point_checked = meshio.read(output / point_file)
        assert np.array_equal(checked.cells_dict[kind], mapping[cells[active]])
        assert np.array_equal(checked.cell_data['material_id'][0], np.asarray(state.get('element_material_ids') or mesh['material_ids'])[active])
        assert np.array_equal(point_checked.point_data['material_id'], point_fields['material_id'])
        assert np.allclose(checked.point_data['displacement_m'], (u - reference)[used])
        for k, v in raw.items():
            assert np.allclose(point_checked.point_data[k], v[qp], equal_nan=True)
        assert np.allclose(point_checked.point_data['displacement_m'], point_u)
        for collection, file in zip(collections, (filename, point_file)):
            ET.SubElement(collection, 'DataSet', timestep=str(frame_index), group=frame['name'], part='0', file=file)
        row = {'stage_index': frame_index, 'stage': frame['name'], 'active_elements': int(active.sum()), 'active_material_points': len(qp), 'relative_force_residual': frame['relative_force_residual']}
        if 'hydrostatic_water_level_m' in state:
            row['hydrostatic_water_level_m'] = h
        for name, m in frame['monitors'].items():
            row[name + '_u_z_m'] = m['mean_displacement_z_m']
            row[name + '_increment_z_m'] = m['stage_increment_z_m']
            row[name + '_minimum_u_z_m'] = m['minimum_displacement_z_m']
            row[name + '_maximum_u_z_m'] = m['maximum_displacement_z_m']
        rows.append(row)
    for root, name in zip(roots, ('stages.pvd', 'integration-points.pvd')):
        ET.ElementTree(root).write(output / name, encoding='utf-8', xml_declaration=True)
    for name, data in [('stages.csv', rows), ('attempts.csv', result['attempts'])]:
        with (output / name).open('x', encoding='utf-8-sig', newline='') as f:
            fields = list(dict.fromkeys((k for row in data for k in row))) if data else ['stage', 'accepted']
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(data)
    if result.get('geostatic') is not None:
        data = result['geostatic']['attempts']
        with (output / 'geostatic-attempts.csv').open('x', encoding='utf-8-sig', newline='') as f:
            fields = list(dict.fromkeys((k for row in data for k in row))) if data else ['stage', 'accepted']
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(data)
        g = result['geostatic']
        summary = {k: g[k] for k in ('complete', 'progress', 'maximum_absolute_engineering_strain_from_seed', 'final_relative_physical_force_residual')}
        summary['accepted_increments'] = sum((a['accepted'] for a in data))
        summary['rejected_increments'] = sum((not a['accepted'] for a in data))
        (output / 'geostatic-summary.json').write_text(json.dumps(summary, indent=2))
    if result['case'].get('joint_system'):
        export_joint_fields(result, output)
    if result['case'].get('beam_system'):
        from beam_io import export
        export(result, output)
        if result['case']['beam_system'].get('plates'):
            from plate_io import export as export_plates
            export_plates(result, output)
    if result['case'].get('anchor_system'):
        from anchor_io import export
        export(result, output)
    if result.get('flow_audits'):
        from export_darcy import export
        for audit in result['flow_audits']:
            hydraulic_mesh = audit['joint_flow']['hydraulic_mesh'] if audit.get('joint_flow') else mesh
            export(hydraulic_mesh, audit['solution'], output / f"flow-{audit['stage_index']:03}")
            if audit.get('joint_flow'):
                with (output / f"joint-flow-map-{audit['stage_index']:03}.json").open('x', encoding='utf-8') as stream:
                    json.dump({k: audit['joint_flow'][k] for k in ('original_mesh_sha256', 'original_case_sha256', 'original_node_to_hydraulic_node', 'original_field')}, stream, indent=2)
    assert fingerprints == {str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (mesh_file, result_file)}, 'export inputs changed'
    record = {'independent_field_accuracy_checked_by_exporter': False, 'schema': 'soil-fem-joint-birth-export/1', 'result_sha256': hashlib.sha256(Path(result_file).read_bytes()).hexdigest(), 'mesh_sha256': hashlib.sha256(Path(mesh_file).read_bytes()).hexdigest(), 'vtu_roundtrip': 'passed', 'integration_point_roundtrip': 'passed', 'cell_field_meaning': 'volume_average', 'frame_count': len(rows)}
    if result['case'].get('joint_system'):
        record['joint_point_roundtrip'] = 'passed'
    if (result['case'].get('beam_system') or {}).get('plates'):
        record['plate_sample_roundtrip'] = 'passed'
    (output / 'export-verification.json').write_text(json.dumps(record, indent=2))
    print(f'Exported {len(rows)} {kind} frames and raw integration-point frames: {output}')

def export_joint_fields(result, output):
    layout = result['joint_integration_points']
    xyz = np.asarray([p['position_m'] for p in layout])
    weights = np.asarray([p['weight_m2'] for p in layout])
    root = ET.Element('VTKFile', type='Collection', version='0.1', byte_order='LittleEndian')
    collection = ET.SubElement(root, 'Collection')
    rows = []
    for fi, frame in enumerate(result['frames']):
        state = frame['state']
        points = state['joints']['points']
        report = frame['joint_results']
        u = np.asarray(state['displacement_m'])
        reference = np.asarray(state.get('display_reference_displacement_m') or result['reference_displacement_m'])
        fields = {}
        for key in ('normal', 'jump_m', 'traction_kpa', 'plastic_slip_m', 'accumulated_plastic_slip_m', 'plastic_dissipation_kpa_m', 'numerical_dissipation_kpa_m', 'endpoint_work_kpa_m'):
            fields[key] = np.asarray([s[key] for s in points])
        for key in ('active', 'gap_m', 'opening_m', 'normal_pressure_kpa', 'shear_traction_kpa', 'elastic_energy_kpa_m', 'contact_reference_offset_m'):
            fields[key] = np.asarray([p[key] for p in report['points']], dtype=int if key == 'active' else float)
        if report.get('water') is not None:
            for key in ('physical_pore_pressure_kpa', 'effective_contact_pressure_kpa', 'total_normal_pressure_kpa', 'friction_limit_kpa'):
                fields[key] = np.asarray([p['hydraulic'][key] for p in report['points']])
            fields['pressure_source'] = np.asarray([{'soil': 0, 'continuous': 1, 'prescribed': 2}[p['hydraulic']['source']] for p in report['points']])
        lifecycle = state['joints'].get('lifecycle', [])
        if lifecycle:
            by_id = {v['joint_id']: v for v in lifecycle}
            life = [by_id[p['joint_id']] for p in layout]
            fields['joint_generation'] = np.asarray([v['generation'] for v in life], dtype=np.uint64)
            fields['joint_prebirth_tied'] = np.asarray([int(v['tied'] and state['active_elements'][p['owner_a']] and state['active_elements'][p['owner_b']]) for v, p in zip(life, layout)])
            fields['birth_policy_code'] = np.asarray([{None: -1, 'transfer_force': 0, 'stress_free': 1}[v['policy']] for v in life])
            fields['has_birth_reference'] = np.asarray([int(p.get('birth') is not None) for p in points])
            for field, key in [('birth_jump_m', 'jump_m'), ('birth_initial_traction_kpa', 'initial_traction_kpa'), ('birth_elastic_offset_m', 'elastic_offset_m')]:
                fields[field] = np.asarray([(p.get('birth') or {}).get(key, [0.0, 0.0, 0.0]) for p in points])
            fields['birth_initial_energy_kpa_m'] = np.asarray([(p.get('birth') or {}).get('initial_energy_kpa_m', 0.0) for p in points])
            archive = output / f'joint-generations-{fi:03}.json'
            with archive.open('x', encoding='utf-8') as stream:
                json.dump(lifecycle, stream, ensure_ascii=False, indent=2)
            assert json.loads(archive.read_text(encoding='utf-8')) == lifecycle
        fields['status'] = np.asarray([{'open': 0, 'stick': 1, 'slip': 2}[s['status']] for s in points])
        fields['integration_point_id'] = np.arange(len(points))
        fields['integration_weight_m2'] = weights
        fields['owner_a'] = np.asarray([p['owner_a'] for p in layout])
        fields['owner_b'] = np.asarray([p['owner_b'] for p in layout])
        for side in ('a', 'b'):
            fields['side_' + side + '_total_displacement_m'] = np.asarray([np.asarray(p['shape']) @ u[p['nodes_' + side]] for p in layout])
        fields['displacement_m'] = np.asarray([np.asarray(p['shape']) @ ((u - reference)[p['nodes_a']] + (u - reference)[p['nodes_b']]) / 2 for p in layout])
        name = f'joints-{fi:03}.vtu'
        meshio.write_points_cells(output / name, xyz, [('vertex', np.arange(len(points))[:, None])], point_data=fields)
        checked = meshio.read(output / name)
        assert np.array_equal(checked.points, xyz)
        for key, value in fields.items():
            assert np.array_equal(checked.point_data[key], value), key
        ET.SubElement(collection, 'DataSet', timestep=str(fi), group=frame['name'], part='0', file=name)
        row = {'stage': frame['name'], 'active_points': int(fields['active'].sum()), 'open_points': int(((fields['status'] == 0) & (fields['active'] == 1)).sum())}
        row.update({k: report[k] for k in ('elastic_energy_kn_m', 'plastic_dissipation_kn_m', 'numerical_dissipation_kn_m', 'endpoint_work_kn_m')})
        rows.append(row)
        if lifecycle:
            row.update(maximum_generation=int(fields['joint_generation'].max()), prebirth_tied_points=int(fields['joint_prebirth_tied'].sum()))
        if report.get('water') is not None:
            row['maximum_physical_pore_pressure_kpa'] = float(fields['physical_pore_pressure_kpa'].max())
    ET.ElementTree(root).write(output / 'joints.pvd', encoding='utf-8', xml_declaration=True)
    with (output / 'joints.csv').open('x', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if any((f['joint_results'].get('water') is not None for f in result['frames'])):
        with (output / 'joint-water-fields.json').open('x', encoding='utf-8') as stream:
            json.dump(dict(pressure_source_codes={'0': 'soil', '1': 'continuous', '2': 'prescribed'}, material_stress_reference_codes={'0': 'soil effective stress', '1': 'impermeable structure total/material stress'}, scope='Effective contact spring stress and physical pressure are separate. Total normal contact pressure equals their sum; friction uses effective contact pressure. Water force integration is stored separately from GP display samples.'), stream, indent=2)
    if any((f['state']['joints'].get('lifecycle') for f in result['frames'])):
        with (output / 'joint-birth-fields.json').open('x', encoding='utf-8') as stream:
            json.dump(dict(birth_policy_codes={'-1': 'legacy active contact', '0': 'transfer_force', '1': 'stress_free'}, birth_reference_validity='has_birth_reference', archive_pattern='joint-generations-NNN.json', scope='Exact data roundtrip only; mechanical BIRTH acceptance uses its dedicated independent reference.'), stream, indent=2)
