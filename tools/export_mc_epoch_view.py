"""Dedicated MC continuation view, anchored to its original Native workflow.

No construction Run is synthesized. The original Run and the exact accepted
snapshot remain separate authoritative sources. The desktop reader must support
soil-fem-solidcraft-mc-epoch-view/1 explicitly.
"""
import argparse
import json
from pathlib import Path, PurePosixPath
import shutil

import numpy as np
from export_mc_epoch import accepted, fields, read_bytes, decode, sha, digest, write_json
from export_solidcraft_view import boundary, clean, MAX_FILE, MAX_TOTAL


def source(root, name):
    parts = PurePosixPath(name)
    if (not name or '\\' in name or ':' in name or parts.is_absolute() or
            any(p in ('', '.', '..') for p in name.split('/'))):
        raise ValueError('unsafe Native source path')
    path = (root/name).resolve()
    if not path.is_relative_to(root):
        raise ValueError('Native source leaves its workflow')
    return path


def export(job, native_origin, output):
    job, native_origin, output = (Path(q).resolve() for q in (job,native_origin,output))
    if output.exists() or any(output.is_relative_to(q) or q.is_relative_to(output) for q in (job,native_origin)):
        raise ValueError('MC Native export requires a new separate folder')
    mesh, origin, checkpoint, verified, snapshot_data, job_manifest = accepted(job)
    native_manifest_bytes = read_bytes(native_origin/'view/view-manifest.json',256*1024)
    native = decode(native_manifest_bytes)
    if native.get('schema') != 'soil-fem-solidcraft-view/1' or native.get('complete') is not True:
        raise ValueError('MC origin must be a completed ordinary continuum Native workflow')
    audit = decode(read_bytes(native_origin/'view-verification.json',256*1024))
    if (audit.get('passed') is not True or audit.get('input_sha256') != native['input_sha256'] or
            audit.get('view_manifest_sha256') != digest(native_manifest_bytes)):
        raise ValueError('Native origin view audit/identity mismatch')
    immutable = job_manifest['binding']['files_sha256']
    for original, frozen in (('analysis/ground/mesh.json','inputs/mesh.json'),
                             ('analysis/stage-results.json','inputs/origin.json'),
                             ('analysis/solver/soil-fem-3d.exe','solver/soil-fem-3d.exe')):
        if native['job_sources_sha256'].get(original) != immutable[frozen]:
            raise ValueError('MC origin mesh/Run/executable differs from the Native workflow')
    if native['solver_binary_sha256'] != immutable['solver/soil-fem-3d.exe']:
        raise ValueError('MC origin solver identity differs')
    copied = dict(native['job_sources_sha256'])
    copied.update({'view/'+n:h for n,h in native['files_sha256'].items()})
    for n in ('view/view-manifest.json','view-verification.json'):
        copied[n] = sha(source(native_origin,n))
    if any(sha(source(native_origin,n)) != h for n,h in copied.items()):
        raise ValueError('Native origin source/packet checksum mismatch')
    geometry = decode(read_bytes(source(native_origin,'view/'+native['geometry']),MAX_FILE))
    cells = np.asarray(mesh.get('tetrahedra10') or mesh['tetrahedra'])
    if (geometry['nodes_m'] != mesh['nodes'] or geometry['cells'] != cells.tolist() or
            geometry['cell_kind'] != ('tetra10' if cells.shape[1]==10 else 'tetra4')):
        raise ValueError('Native geometry differs from the fixed MC mesh')
    layout = origin['integration_points']
    owners = np.asarray([r['element_index'] for r in layout])
    current = checkpoint['payload']['material_checkpoint']['payload']['committed']
    continuation = checkpoint['payload']['continuation']
    raw_parameters = {'young_kpa','poisson','cohesion_kpa','friction_deg','dilation_deg','tension_cutoff_kpa','has_tension_cutoff'}
    prepared = []
    for index,name in enumerate(('Unreduced construction endpoint','Accepted material continuation')):
        state = current if index else origin['last_committed']
        _,_,_,nodal,averages,positions,gp = fields(mesh,origin,checkpoint,bool(index))
        active = np.asarray(state['active_elements'],dtype=bool)
        ids = np.flatnonzero(active);used=nodal['original_node_id'];qp=gp['integration_point_id']
        cell_fields = {n.removeprefix('volume_average_'):v[0] for n,v in averages.items()
                       if n.removeprefix('volume_average_') not in raw_parameters and n!='original_element_id'}
        cell_fields['material_model'] = np.ones(len(ids),dtype=int)
        node_fields = {n:v for n,v in nodal.items() if n!='original_node_id'}
        point_fields = {n:v for n,v in gp.items() if n not in raw_parameters | {
            'integration_point_id','original_element_id','integration_weight_m3',
            'epoch_displacement_increment_m','total_displacement_from_seed_m'}}
        metrics = (continuation['accepted'][-1] if index and continuation['accepted'] else origin['frames'][-1])
        frame = clean(dict(schema='soil-fem-view-frame/4',index=index,name=name,
            active_element_ids=ids,active_node_ids=used,
            boundary_faces=boundary(cells,active,np.asarray(mesh['nodes'])),node_fields=node_fields,
            cell_fields=cell_fields,monitors={},relative_force_residual=metrics['relative_force_residual'],
            maximum_normalized_yield=metrics['maximum_normalized_yield'],
            field_semantics=dict(source_family='mc_parameter_continuation',
                increment='accepted endpoint minus the original unreduced construction endpoint')))
        points = clean(dict(schema='soil-fem-view-integration-points/4',frame_index=index,
            integration_point_ids=qp,element_ids=owners[qp],positions_m=positions,
            weights_m3=gp['integration_weight_m3'],fields=point_fields))
        prepared.append((frame,points))
    encoded = [json.dumps(v,ensure_ascii=False,allow_nan=False).encode('utf-8') for pair in prepared for v in pair]
    if any(len(v)>MAX_FILE for v in encoded) or sum(map(len,encoded))>MAX_TOTAL:
        raise ValueError('MC view exceeds unchanged display budgets')
    output.mkdir(parents=True,exist_ok=False);view=output/'view';view.mkdir()
    source_hashes={}
    def copy(src,name):
        target=output/name;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(src,target);source_hashes[name]=sha(target)
        if source_hashes[name]!=sha(src):raise ValueError('MC source changed while copying')
    for n in sorted(copied):copy(source(native_origin,n),'origin/'+n)
    for n,target in (('inputs/mesh.json','epoch/mesh.json'),('inputs/origin.json','epoch/origin.json'),('job.json','epoch/job.json')):
        copy(job/n,target)
    (output/'epoch/accepted-snapshot.json').write_bytes(snapshot_data)
    write_json(output/'epoch/verified.json',verified)
    for n in ('epoch/accepted-snapshot.json','epoch/verified.json'):source_hashes[n]=sha(output/n)
    shutil.copyfile(source(native_origin,'view/'+native['geometry']),view/'geometry.json')
    files={'geometry.json':sha(view/'geometry.json')};entries=[]
    for i,(frame,points) in enumerate(prepared):
        for stem,value in (('frame',frame),('points',points)):
            name=f'{stem}-{i:03}.json';write_json(view/name,value);files[name]=sha(view/name)
        entries.append(dict(index=i,name=frame['name'],frame=f'frame-{i:03}.json',integration_points=f'points-{i:03}.json'))
    metadata=dict(origin_manifest='origin/view/view-manifest.json',snapshot='epoch/accepted-snapshot.json',
        origin='epoch/origin.json',mesh='epoch/mesh.json',job='epoch/job.json',verification='epoch/verified.json',
        parameter_progress=verified['progress'],path_complete=verified['complete'],stop_reason=verified['stop_reason'],
        stability_classification='not_evaluated')
    manifest=dict(schema='soil-fem-solidcraft-mc-epoch-view/1',complete=True,
        input_sha256=native['input_sha256'],cad_sha256=native['cad_sha256'],
        exporter_binary_sha256=native['exporter_binary_sha256'],solver_binary_sha256=native['solver_binary_sha256'],
        geometry='geometry.json',frames=entries,files_sha256=files,job_sources_sha256=source_hashes,
        semantics=dict(length_unit='m',stress_unit='kPa',stress_sign='compression positive',
            stress_components=['xx','yy','zz','xy','yz','xz'],source_result_schema='soil-fem-mc-epoch-snapshot/1',
            path='two verified accepted endpoints; incomplete is not collapse',
            epoch_displacement_increment_m='current minus unreduced construction endpoint; z up'),
        omissions=['water/structures unsupported by the MC continuation solver','intermediate material histories'],
        scope='Dedicated dry continuum MC continuation; not an inferred safety factor',epoch=metadata)
    if any(sha(source(native_origin,n))!=h for n,h in copied.items()):
        raise ValueError('Native origin changed while exporting MC packet')
    for n in ('inputs/mesh.json','inputs/origin.json'):
        if sha(job/n)!=immutable[n]:raise ValueError('MC input changed while exporting')
    write_json(view/'view-manifest.json',manifest)
    write_json(output/'view-verification.json',dict(schema='soil-fem-mc-epoch-view-export/1',passed=True,
        input_sha256=native['input_sha256'],view_manifest_sha256=sha(view/'view-manifest.json'),
        mechanical_revalidation=verified,source_snapshot_exact=True,ordinary_run_not_synthesized=True,
        native_reader_acceptance_pending=True,p6_complete=False))
    return manifest


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('job',type=Path)
    parser.add_argument('native_origin',type=Path);parser.add_argument('output',type=Path)
    a=parser.parse_args();result=export(a.job,a.native_origin,a.output)
    print(json.dumps(dict(schema=result['schema'],frames=len(result['frames']),**result['epoch'])))
