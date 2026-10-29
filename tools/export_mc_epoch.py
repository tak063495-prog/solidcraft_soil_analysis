"""Export a deeply verified, accepted MC continuation snapshot without changing it.

The dedicated packet is not an ordinary construction Run or a stability result.
It includes endpoint VTUs, raw GP VTUs and the exact accepted checkpoint.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

import meshio
import numpy as np
from high_order_io import state_fields

LIMIT = 1024**3
FILES = {'inputs/mesh.json', 'inputs/origin.json', 'inputs/specification.json',
         'solver/soil-fem-3d.exe'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_bytes(path, maximum=LIMIT):
    with path.open('rb') as stream:
        data = stream.read(maximum + 1)
    if len(data) > maximum:
        raise ValueError('MC export input exceeds size limit: ' + str(path))
    return data


def decode(data):
    def reject(value):
        raise ValueError('non-finite JSON number: ' + value)
    def real(value):
        result = float(value)
        if not math.isfinite(result):
            reject(value)
        return result
    return json.loads(data, parse_constant=reject, parse_float=real)


def accepted(job):
    """Bind in-memory bytes to the same immutable files and verified snapshot.

    A running writer may advance either alternating slot. Retry at most three
    times; never present a mixture of generations or trust report/status liveness.
    """
    manifest_data = read_bytes(job / 'job.json', 256 * 1024)
    manifest = decode(manifest_data)
    if manifest['schema'] != 'soil-fem-mc-epoch-job/1':
        raise ValueError('expected a dedicated MC continuation job')
    hashes = manifest['binding']['files_sha256']
    expected = FILES | ({'inputs/checkpoint.json'} if manifest['binding']['parent'] else set())
    if set(hashes) != expected:
        raise ValueError('unexpected MC immutable inventory')
    immutable = {n: read_bytes(job / n) for n in hashes if n != 'solver/soil-fem-3d.exe'}
    if any(digest(data) != hashes[n] for n, data in immutable.items()):
        raise ValueError('MC immutable input checksum mismatch')
    solver = job / 'solver/soil-fem-3d.exe'
    if sha(solver) != hashes['solver/soil-fem-3d.exe']:
        raise ValueError('MC frozen executable checksum mismatch')
    for _ in range(3):
        pointer_data = read_bytes(job / 'recovery.json', 8192)
        pointer = decode(pointer_data)
        generation = pointer['generation']
        if (type(generation) is not int or generation < 1 or
                pointer['schema'] != 'soil-fem-mc-epoch-recovery/1' or
                pointer['job_id'] != manifest['job_id'] or
                pointer['checkpoint'] != f'recovery-checkpoint-{generation % 2}.json'):
            raise ValueError('MC recovery pointer identity/path mismatch')
        snapshot_data = read_bytes(job / pointer['checkpoint'])
        if digest(snapshot_data) != pointer['checkpoint_sha256']:
            if read_bytes(job / 'recovery.json', 8192) != pointer_data:
                continue
            raise ValueError('MC accepted snapshot checksum mismatch')
        checked = subprocess.run([str(solver), 'verify-mc-epoch-job', str(job)],
            capture_output=True, timeout=1800,
            creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
        if (read_bytes(job / 'recovery.json', 8192) != pointer_data or
                sha(job / pointer['checkpoint']) != digest(snapshot_data)):
            continue
        if checked.returncode:
            raise ValueError('MC mechanical checkpoint verification failed: ' +
                             checked.stderr.decode('utf-8', errors='replace')[-4000:])
        info = decode(checked.stdout)
        snapshot = decode(snapshot_data)
        checkpoint = snapshot['checkpoint']
        if (info['job_id'] != manifest['job_id'] or info['generation'] != generation or
                info['checkpoint_sha256'] != digest(snapshot_data) or
                info['payload_sha256'] != checkpoint['payload_sha256'] or
                pointer['payload_sha256'] != checkpoint['payload_sha256'] or
                snapshot['job_id'] != manifest['job_id'] or snapshot['generation'] != generation):
            raise ValueError('MC verifier/snapshot generation binding mismatch')
        if (read_bytes(job / 'job.json', 256 * 1024) != manifest_data or
                any(sha(job / n) != h for n, h in hashes.items())):
            raise ValueError('MC immutable source changed during verification')
        return (decode(immutable['inputs/mesh.json']), decode(immutable['inputs/origin.json']),
                checkpoint, info, snapshot_data, manifest)
    raise ValueError('MC recovery advanced during all three reads; retry export')


def fields(mesh, origin, checkpoint, current):
    payload = checkpoint['payload']['material_checkpoint']['payload']
    state = payload['committed'] if current else origin['last_committed']
    xyz = np.asarray(mesh['nodes'], dtype=float)
    cells = np.asarray(mesh.get('tetrahedra10') or mesh['tetrahedra'], dtype=int)
    kind = 'tetra10' if cells.shape[1] == 10 else 'tetra'
    active = np.asarray(state['active_elements'], dtype=bool)
    ids = np.flatnonzero(active)
    if not len(ids):
        raise ValueError('MC snapshot has no active continuum cells to display')
    used = np.unique(cells[ids])
    remap = np.full(len(xyz), -1, dtype=int)
    remap[used] = np.arange(len(used))
    u = np.asarray(state['displacement_m'], dtype=float)
    base = np.asarray(origin['last_committed']['displacement_m'], dtype=float)
    reference = np.asarray(state.get('display_reference_displacement_m') or
                           origin.get('reference_displacement_m', np.zeros_like(xyz)), dtype=float)
    nodal = dict(original_node_id=used, total_displacement_from_seed_m=u[used],
                 displacement_m=(u-reference)[used], display_reference_displacement_m=reference[used],
                 epoch_displacement_increment_m=(u-base)[used])
    layout = origin['integration_points']
    owners = np.asarray([p['element_index'] for p in layout], dtype=int)
    qp = np.flatnonzero(active[owners])
    weights = np.asarray([p['weight_m3'] for p in layout], dtype=float)
    grouped = [[] for _ in cells]
    for qi, ei in enumerate(owners):
        grouped[ei].append(qi)
    groups = [np.asarray(grouped[ei], dtype=int) for ei in ids]
    raw = state_fields(state['material_states'])
    materials = np.asarray(state.get('element_material_ids') or mesh['material_ids'])
    # MATSET uses material IDs, not table offsets. Removed elements are not
    # exported; their old histories stay verbatim in accepted-snapshot.json.
    params = payload['parameters'] if current else [m['mohr_coulomb'] for m in origin['ground']['materials']]
    by_id = dict(zip([m['id'] for m in origin['ground']['materials']], params))
    material_fields = {name: np.asarray([by_id[int(materials[ei])][name] for ei in owners])
                       for name in ('young_kpa', 'poisson', 'cohesion_kpa', 'friction_deg', 'dilation_deg')}
    material_fields['has_tension_cutoff'] = np.asarray([
        int(by_id[int(materials[ei])]['tension_cutoff_kpa'] is not None) for ei in owners])
    # NaN means unavailable on MC; the explicit mask distinguishes it from zero.
    material_fields['tension_cutoff_kpa'] = np.asarray([
        by_id[int(materials[ei])]['tension_cutoff_kpa'] if
        by_id[int(materials[ei])]['tension_cutoff_kpa'] is not None else np.nan for ei in owners])
    raw.update(material_fields)
    cell_fields = {('volume_average_' + n): [np.asarray([
        np.average(values[g], axis=0, weights=weights[g]) for g in groups])] for n, values in raw.items()}
    cell_fields.update(original_element_id=[ids], material_id=[materials[ids]],
                       region_id=[np.asarray(mesh['region_ids'])[ids]])
    gp_u = np.asarray([np.asarray(layout[i]['shape_functions']) @ u[cells[owners[i]]] for i in qp])
    gp_du = np.asarray([np.asarray(layout[i]['shape_functions']) @ (u-base)[cells[owners[i]]] for i in qp])
    gp_fields = {n: values[qp] for n, values in raw.items()}
    gp_fields.update(integration_point_id=qp, original_element_id=owners[qp],
        material_id=materials[owners[qp]], integration_weight_m3=weights[qp],
        total_displacement_from_seed_m=gp_u, epoch_displacement_increment_m=gp_du)
    return (xyz[used], kind, remap[cells[ids]], nodal, cell_fields,
            np.asarray([layout[i]['position_m'] for i in qp]), gp_fields)


def write_json(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, allow_nan=False, indent=2)


def export(job, output):
    job, output = Path(job).resolve(), Path(output).resolve()
    if output.exists() or output.is_relative_to(job) or job.is_relative_to(output):
        raise ValueError('MC export requires a new folder outside the source job')
    mesh, origin, checkpoint, verified, snapshot_data, manifest = accepted(job)
    endpoints = {name: fields(mesh, origin, checkpoint, name == 'current') for name in ('origin', 'current')}
    output.mkdir(parents=True, exist_ok=False)
    for name, (xyz, kind, cells, nodal, cell, gp_xyz, gp) in endpoints.items():
        meshio.write_points_cells(output / (name + '.vtu'), xyz, [(kind, cells)],
                                  point_data=nodal, cell_data=cell)
        meshio.write_points_cells(output / (name + '-points.vtu'), gp_xyz,
            [('vertex', np.arange(len(gp_xyz))[:, None])], point_data=gp)
        restored = meshio.read(output / (name + '.vtu'))
        restored_gp = meshio.read(output / (name + '-points.vtu'))
        if (not np.array_equal(restored.points, xyz) or
                not np.array_equal(restored.cells_dict[kind], cells) or
                not np.array_equal(restored_gp.points, gp_xyz)):
            raise ValueError('MC export geometry roundtrip failed')
        for expected, actual in ((nodal, restored.point_data), (gp, restored_gp.point_data)):
            if set(expected) != set(actual) or any(not np.array_equal(v, actual[n], equal_nan=True) for n,v in expected.items()):
                raise ValueError('MC export point-field roundtrip failed')
        if set(cell) != set(restored.cell_data) or any(
                not np.array_equal(v[0], restored.cell_data[n][0], equal_nan=True) for n,v in cell.items()):
            raise ValueError('MC export cell-field roundtrip failed')
    with (output / 'accepted-snapshot.json').open('xb') as stream:
        stream.write(snapshot_data)
    payload = checkpoint['payload']['material_checkpoint']['payload']
    packet = dict(schema='soil-fem-mc-epoch-export/1', job_id=manifest['job_id'],
        verified_snapshot=verified, mechanical_revalidation='passed', vtu_roundtrip='passed',
        complete=verified['complete'], stability_classification='not_evaluated',
        parameter_progress=verified['progress'], stop_reason=verified['stop_reason'],
        source_files_sha256=manifest['binding']['files_sha256'],
        current_parameters=payload['parameters'],
        origin_parameters=[m['mohr_coulomb'] for m in origin['ground']['materials']],
        material_table_ids=[m['id'] for m in origin['ground']['materials']],
        files_sha256={p.name:sha(p) for p in output.iterdir() if p.is_file()},
        field_semantics=dict(coordinates='undeformed metres; z positive upward',
            total_displacement_from_seed_m='cumulative constitutive displacement; no display reset',
            displacement_m='total minus the retained display reference',
            epoch_displacement_increment_m='current minus the unreduced construction endpoint; negative z settlement, positive z heave',
            stress='compression-positive effective stress at raw GP; volume averages identified separately',
            inactive='removed cells omitted from VTU; every frozen history retained in accepted-snapshot.json',
            parameters='current parameters on active cells only; null Rankine cutoff represented by NaN with has_tension_cutoff=0',
            path='two endpoint states only; intermediate accepted epoch histories are not reconstructed',
            stability='lambda is a material path coordinate; incomplete is not a physical collapse bound'),
        native_solidcraft_reader_supported=False)
    write_json(output / 'export-verification.json', packet)
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('job', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = export(args.job, args.output)
    print(json.dumps({k:result[k] for k in ('schema','complete','parameter_progress','vtu_roundtrip')}))
