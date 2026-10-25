"""Primary cell fields and raw RT0 quadrature values, no nodal smoothing."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import meshio
import numpy as np
import fem_reference as fem
from richards_io import geometry,fields


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def export(mesh_path,result_path,solver,output):
    result=json.loads(result_path.read_text(encoding='utf-8-sig'));mesh=json.loads(mesh_path.read_text(encoding='utf-8-sig'))
    assert result['schema']=='soil-fem-richards-result/1' and result['complete'],'only a complete hydraulic result can be exported'
    assert sha(solver)==result['solver_binary_sha256'],'export requires the result executable'
    output.mkdir(parents=True,exist_ok=False)
    embedded=output/'identity-case.json';embedded.write_text(json.dumps(result['case']),encoding='utf-8')
    command=[str(solver),'richards-identity',str(mesh_path),str(embedded),str(output/'identity.json')]
    log=subprocess.run(command,capture_output=True,text=True,encoding="utf-8",errors="replace")
    (output/'identity.log').write_text(log.stdout+log.stderr,encoding='utf-8');assert log.returncode==0
    identity=json.loads((output/'identity.json').read_text(encoding='utf-8-sig'));assert identity['mesh_sha256']==result['mesh_sha256'] and identity['case_sha256']==result['case_sha256']
    xyz=np.asarray(mesh['nodes']);cells=fem.cells(mesh);active=np.flatnonzero(result['active_elements']);used=np.unique(cells[active])
    mapping=np.full(len(xyz),-1);mapping[used]=np.arange(len(used));kind='tetra10' if cells.shape[1]==10 else 'tetra'
    geometries={e:geometry(mesh,e) for e in active};root=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian');collection=ET.SubElement(root,'Collection')
    raw_root=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian');raw_collection=ET.SubElement(raw_root,'Collection');rows=[]
    maximum_flux_error=maximum_storage_error=0.
    for i,frame in enumerate(result['frames']):
        blocks,audit=fields(mesh,result,frame,geometries)
        maximum_flux_error=max(maximum_flux_error,audit['maximum_independent_local_flux_error_m3_s']);maximum_storage_error=max(maximum_storage_error,audit['maximum_independent_storage_error_m3'])
        data={'original_element_id':[active],'material_id':[np.asarray(result['material_ids'])[active]],
            'cell_total_head_m':[np.asarray([result['case']['head_datum_m']+frame['state']['relative_cell_head_m'][e] for e in active])],
            'stored_water_offset_m3':[np.asarray([frame['state']['cell_storage_offset_m3'][e] for e in active])]}
        names={'pressure_head':'pressure_head_m','water_content':'water_content','effective_saturation':'effective_saturation','water_saturation':'water_saturation','relative_conductivity':'relative_conductivity','flux':'darcy_flux_m_s'}
        for name,label in names.items():
            data['volume_average_'+label]=[np.asarray([np.average(b[name],weights=b['weight'],axis=0) for b in blocks])]
        filename=f'hydraulic-{i:04}.vtu';meshio.write_points_cells(output/filename,xyz[used],[(kind,mapping[cells[active]])],point_data={'original_node_id':used},cell_data=data)
        checked=meshio.read(output/filename);assert np.array_equal(checked.cells_dict[kind],mapping[cells[active]])
        for name,v in data.items():assert np.array_equal(checked.cell_data[name][0],v[0])
        raw={label:np.concatenate([b[name] for b in blocks]) for name,label in {'head':'total_head_m','weight':'integration_weight_m3',**names}.items()}
        raw['original_element_id']=np.concatenate([np.full(len(b['weight']),b['element']) for b in blocks])
        points=np.vstack([b['position'] for b in blocks]);raw_name=f'hydraulic-points-{i:04}.vtu'
        meshio.write_points_cells(output/raw_name,points,[('vertex',np.arange(len(points))[:,None])],point_data=raw)
        checked_raw=meshio.read(output/raw_name)
        for name,value in raw.items():assert np.array_equal(checked_raw.point_data[name],value)
        if result['case']['analysis']['kind']!='steady' or i==len(result['frames'])-1:
            ET.SubElement(collection,'DataSet',timestep=str(frame['state']['time_s']),group='',part='0',file=filename)
            ET.SubElement(raw_collection,'DataSet',timestep=str(frame['state']['time_s']),group='',part='0',file=raw_name)
        rows.append(dict(frame=i,name=frame['name'],time_s=frame['state']['time_s'],**{k:v for k,v in frame['balance'].items() if k!='atmospheric_faces'}))
    ET.ElementTree(root).write(output/'hydraulics.pvd',encoding='utf-8',xml_declaration=True);ET.ElementTree(raw_root).write(output/'hydraulic-points.pvd',encoding='utf-8',xml_declaration=True)
    with (output/'water-balance.csv').open('x',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    with (output/'export-verification.json').open('x',encoding='utf-8') as stream:json.dump(dict(schema='soil-fem-richards-export/1',passed=True,frames=len(rows),cell_and_raw_point_roundtrip='passed',
        maximum_independent_local_flux_error_m3_s=maximum_flux_error,maximum_independent_storage_error_m3=maximum_storage_error,
        mesh_file_sha256=sha(mesh_path),result_file_sha256=sha(result_path),solver_sha256=sha(solver),
        scope='RT0 primary face/cell fields and reconstructed original quadrature values. No nodal pressure/saturation smoothing. Constitutive reassembly is an audit, not an independent nonlinear PDE solve.'),stream,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('result',type=Path);p.add_argument('output',type=Path);p.add_argument('--solver',type=Path,required=True);a=p.parse_args()
    export(a.mesh.resolve(),a.result.resolve(),a.solver.resolve(),a.output.resolve())
