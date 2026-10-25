"""Frame lines, physical nodal kinematics and unsmoothed section actions."""
import csv
import json
import xml.etree.ElementTree as ET
import meshio
import numpy as np


def empty_frame(path):
    # Empty frames must be in the time series, otherwise removed members can
    # remain visible at the previous timestep. meshio cannot read zero cells.
    root=ET.Element('VTKFile',type='UnstructuredGrid',version='0.1',byte_order='LittleEndian')
    piece=ET.SubElement(ET.SubElement(root,'UnstructuredGrid'),'Piece',NumberOfPoints='0',NumberOfCells='0')
    ET.SubElement(ET.SubElement(piece,'Points'),'DataArray',type='Float64',NumberOfComponents='3',format='ascii')
    cells=ET.SubElement(piece,'Cells')
    for name,kind in (('connectivity','Int64'),('offsets','Int64'),('types','UInt8')):
        ET.SubElement(cells,'DataArray',type=kind,Name=name,format='ascii')
    ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
    assert ET.parse(path).find('UnstructuredGrid/Piece').get('NumberOfCells')=='0'


def export(run,output):
    cfg=run['case']['beam_system'];node_ids={n['id']:i for i,n in enumerate(cfg['nodes'])}
    member_ids={m['id']:i for i,m in enumerate(cfg['members'])}
    root=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian');collection=ET.SubElement(root,'Collection')
    node_rows=[];member_rows=[]
    for fi,frame in enumerate(run['frames']):
        report=frame['beam_results'];state=frame['state']['beams']
        assert np.array_equal(report['node_displacement_m'],state['node_displacement_m'])
        assert np.array_equal(report['node_rotation_rad'],state['node_rotation_rad'])
        used=np.flatnonzero(report['node_active']);mapping={n:i for i,n in enumerate(used)}
        filename=f'beams-{fi:03}.vtu';members=report['members']
        assert len(members)==sum(state['active_members'])
        if not members:empty_frame(output/filename)
        else:
            xyz=np.asarray(report['node_positions_m']);u=np.asarray(report['node_displacement_m']);rotation=np.asarray(report['node_rotation_rad'])
            indices=[member_ids[m['id']] for m in members]
            lines=np.array([[mapping[node_ids[n]] for n in cfg['members'][mi]['nodes']] for mi in indices])
            actions=np.asarray([m['local_end_actions'] for m in members]);section=np.asarray([m['section_forces_at_start_mid_end'] for m in members])
            fields={'member_index':np.asarray(indices),'length_m':np.array([m['length_m'] for m in members])}
            for end,start in (('i',0),('j',6)):
                fields[f'end_{end}_force_local_kn']=actions[:,start:start+3]
                fields[f'end_{end}_moment_local_kn_m']=actions[:,start+3:start+6]
            for j,name in enumerate(('start','mid','end')):
                fields[f'section_{name}_force_local_kn']=section[:,j,:3]
                fields[f'section_{name}_moment_local_kn_m']=section[:,j,3:]
            for j,name in enumerate(('x','y','z')):fields[f'local_axis_{name}']=np.asarray([m['local_axes'][j] for m in members])
            points={'node_index':used,'total_displacement_from_cad_m':u[used],'total_rotation_rad':rotation[used]}
            meshio.write_points_cells(output/filename,xyz[used],[('line',lines)],point_data=points,cell_data={k:[v] for k,v in fields.items()})
            check=meshio.read(output/filename)
            assert np.array_equal(check.cells_dict['line'],lines)
            for k,v in points.items():assert np.allclose(check.point_data[k],v)
            for k,v in fields.items():assert np.allclose(check.cell_data[k][0],v)
        ET.SubElement(collection,'DataSet',timestep=str(fi),group=frame['name'],part='0',file=filename)
        for i,node in enumerate(cfg['nodes']):
            row={'stage_index':fi,'stage':frame['name'],'node_id':node['id'],'active':report['node_active'][i]}
            for j,a in enumerate('xyz'):
                row[f'cad_{a}_m']=report['node_positions_m'][i][j]
                row[f'total_u_{a}_m']=state['node_displacement_m'][i][j]
                row[f'total_rotation_{a}_rad']=state['node_rotation_rad'][i][j]
            node_rows.append(row)
        for m in members:
            row={'stage_index':fi,'stage':frame['name'],'member_id':m['id'],'length_m':m['length_m']}
            for j,name in enumerate(('start','mid','end')):
                for a,value in zip(('N_kn','Vy_kn','Vz_kn','T_kn_m','My_kn_m','Mz_kn_m'),m['section_forces_at_start_mid_end'][j]):row[f'{name}_{a}']=value
            member_rows.append(row)
    ET.ElementTree(root).write(output/'beams.pvd',encoding='utf-8',xml_declaration=True)
    for name,rows in (('beam-nodes.csv',node_rows),('beam-members.csv',member_rows)):
        with (output/name).open('x',encoding='utf-8-sig',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['stage_index','member_id'])
            writer.writeheader();writer.writerows(rows)
    metadata={'schema':'soil-fem-beam-export/1','nodes':node_ids,'members':member_ids,
              'geometry':'Undeformed CAD chord; warp with total_displacement_from_cad_m. Physical beam kinematics use the CAD origin, independent of the soil display RESET.',
              'section_actions':'No smoothing. N is tension positive. Local axes and i/j element end actions are recorded separately.',
              'roundtrip':'Nonempty frames read with meshio; zero-cell frames checked as explicit VTK XML.'}
    (output/'beam-fields.json').write_text(json.dumps(metadata,indent=2))
