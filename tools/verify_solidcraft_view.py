"""Read a scene independently against the job's existing VTU exports.

This checks display values/topology and provenance, not a new physical solution
or a native GUI. Synthetic CAD fixtures remain synthetic after this audit.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import meshio
import numpy as np
from export_solidcraft_view import read,read_result
from verify_acceptance_packet import sha


def verify(job,packet,output):
    job=job.resolve();packet=packet.resolve();manifest=read(packet/'view-manifest.json')
    if manifest.get('schema') not in ('soil-fem-solidcraft-view/1','soil-fem-solidcraft-view/2','soil-fem-solidcraft-view/3') or manifest.get('complete') is not True:raise ValueError('complete view packet required')
    structural=manifest['schema']=='soil-fem-solidcraft-view/3'
    projected=manifest['schema']=='soil-fem-solidcraft-view/2' or manifest['semantics'].get('source_result_schema')=='soil-fem-projected-mcc-stage-result/2'
    if projected and manifest['semantics'].get('source_result_schema') not in ('soil-fem-projected-mcc-stage-result/1','soil-fem-projected-mcc-stage-result/2'):raise ValueError('projected source schema required')
    for base,files in ((job,manifest['job_sources_sha256']),(packet,manifest['files_sha256'])):
        for name,digest in files.items():
            path=(base/name).resolve()
            if not path.is_relative_to(base) or sha(path)!=digest:raise ValueError('view packet/source file hash mismatch')
    geometry=read(packet/manifest['geometry']);xyz=np.asarray(geometry['nodes_m']);cells=np.asarray(geometry['cells'])
    structural_source=(read(job/'analysis/ground/mesh.json'),read_result(job/'analysis/stage-results.json')) if structural else None
    results=job/'analysis/results';reports=[]
    for frame in manifest['frames']:
        i=frame['index'];view=read(packet/frame['frame']);points=read(packet/frame['integration_points'])
        version='3' if structural else '2' if projected else '1'
        if view.get('schema')!='soil-fem-view-frame/'+version or points.get('schema')!='soil-fem-view-integration-points/'+version:raise ValueError('view frame family differs from manifest')
        vtk=meshio.read(results/f'stage-{i:03}.vtu');kind='tetra10' if geometry['cell_kind']=='tetra10' else 'tetra'
        node_ids=np.asarray(view['active_node_ids']);element_ids=np.asarray(view['active_element_ids'])
        if not np.array_equal(node_ids,vtk.point_data['original_node_id']):raise ValueError('original node identifiers differ from VTU')
        if not np.array_equal(element_ids,vtk.cell_data['original_element_id'][0]):raise ValueError('original element identifiers differ from VTU')
        if not np.array_equal(xyz[node_ids],vtk.points) or not np.array_equal(cells[element_ids],node_ids[vtk.cells_dict[kind]]):raise ValueError('volume geometry differs from VTU')
        checked_node=[];checked_cell=[]
        for name,values in view['node_fields'].items():
            if name not in vtk.point_data:continue
            if not np.allclose(values,vtk.point_data[name],rtol=1e-12,atol=1e-14):raise ValueError('nodal view field differs from VTU: '+name)
            checked_node.append(name)
        for name,values in view['cell_fields'].items():
            key=name if name in vtk.cell_data else 'volume_average_'+name
            if name=='material_model' and key not in vtk.cell_data and kind=='tetra':
                if 'preconsolidation_pressure_kpa' not in vtk.cell_data or not np.all(np.asarray(values)==0):raise ValueError('legacy MCC model metadata mismatch')
                continue
            if key not in vtk.cell_data:raise ValueError('missing cell reference in VTU: '+name)
            if not np.allclose(np.asarray(values,dtype=float),vtk.cell_data[key][0],rtol=1e-12,atol=1e-12,equal_nan=True):raise ValueError('cell view field differs from VTU: '+name)
            checked_cell.append(name)
        faces=Counter(tuple(sorted(int(cell[a]) for a in local)) for cell in cells[element_ids]
            for local in ((0,1,2),(0,1,3),(0,2,3),(1,2,3)))
        expected={key for key,count in faces.items() if count==1}
        actual=[tuple(sorted(f['nodes'][:3])) for f in view['boundary_faces']]
        if len(set(actual))!=len(actual) or set(actual)!=expected:raise ValueError('display boundary differs from active VTU topology')
        raw_file=results/f'integration-points-{i:03}.vtu';checked_raw=[]
        if raw_file.exists():
            raw=meshio.read(raw_file)
            if not np.array_equal(points['integration_point_ids'],raw.point_data['integration_point_id']):raise ValueError('raw point identifiers differ from VTU')
            if not np.array_equal(points['positions_m'],raw.points):raise ValueError('raw point geometry differs from VTU')
            if not np.array_equal(points['element_ids'],raw.point_data['original_element_id']):raise ValueError('raw point owners differ from VTU')
            if not np.allclose(points['weights_m3'],raw.point_data['integration_weight_m3'],rtol=1e-13,atol=1e-14):raise ValueError('raw point weights differ from VTU')
            for name,values in points['fields'].items():
                if name not in raw.point_data:continue
                if not np.allclose(np.asarray(values,dtype=float),raw.point_data[name],rtol=1e-13,atol=1e-13,equal_nan=True):raise ValueError('raw point view field differs from VTU: '+name)
                checked_raw.append(name)
        else:
            if kind!='tetra':raise ValueError('quadratic raw point VTU is missing')
            expected_points=xyz[cells[element_ids]].mean(axis=1)
            if not np.array_equal(points['element_ids'],element_ids) or not np.array_equal(points['positions_m'],expected_points):raise ValueError('legacy single-point geometry differs from cell centroids')
            tets=xyz[cells[element_ids]];volumes=np.linalg.det(np.transpose(tets[:,1:]-tets[:,:1],(0,2,1)))/6
            if not np.array_equal(points['integration_point_ids'],element_ids) or not np.allclose(points['weights_m3'],volumes,rtol=1e-13,atol=1e-14):raise ValueError('legacy single-point identity/volume mismatch')
            for name,values in points['fields'].items():
                if name not in vtk.cell_data:continue
                if not np.allclose(np.asarray(values,dtype=float),vtk.cell_data[name][0],rtol=1e-13,atol=1e-13,equal_nan=True):raise ValueError('legacy single-point value differs from VTU: '+name)
                checked_raw.append(name)
        if projected:
            primary=view['mixed_primary'];vtk_primary=meshio.read(results/f'mean-pressure-dofs-{i:03}.vtu')
            if not np.array_equal(primary['positions_m'],vtk_primary.points):raise ValueError('primary coordinates differ from VTU')
            for name,key in (('original_node_ids','original_node_id'),('material_ids','material_id'),
                ('mean_effective_stress_compression_positive_kpa','mean_effective_stress_compression_positive_kpa'),
                ('cumulative_material_volume_strain_compression_positive','cumulative_material_volume_strain_compression_positive'),
                ('frozen_compatibility_offset_m3','frozen_compatibility_offset_m3')):
                if not np.array_equal(primary[name],vtk_primary.point_data[key]):raise ValueError('primary view field differs from VTU: '+name)
            if view['field_semantics'].get('equilibrium_stress')!='work_conjugate_effective_stress_compression_positive_kpa':raise ValueError('projected equilibrium stress semantics differ')
        structural_counts=None
        if structural:
            from solidcraft_structure_fields import fields
            from verify_structural_view_fields import verify as verify_structure_vtu
            source_mesh,source_run=structural_source
            if view.get('structures')!=fields(source_mesh,source_run,i):raise ValueError('structural display differs from authoritative source')
            structural_counts=verify_structure_vtu(view['structures'],results,i)
        reports.append(dict(index=i,name=view['name'],active_nodes=len(node_ids),active_elements=len(element_ids),boundary_faces=len(actual),
            raw_points=len(points['element_ids']),checked_node_fields=checked_node,checked_volume_average_fields=checked_cell,checked_raw_fields=checked_raw,
            **({'checked_structural_cells':structural_counts} if structural else {})))
    report=dict(schema='soil-fem-solidcraft-view-verification/1',passed=True,input_sha256=manifest['input_sha256'],
        view_manifest_sha256=sha(packet/'view-manifest.json'),frames=reports,
        scope='External display data compared with independently read existing VTU exports; not a native GUI or new physical reference solution.')
    with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print('Viewer-data/VTU audit passed:',len(reports),'frames')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('job',type=Path);p.add_argument('packet',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();verify(a.job,a.packet,a.output)
