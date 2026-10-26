"""Compare structural viewer transport with separately read VTU output."""
import xml.etree.ElementTree as ET
import meshio
import numpy as np


def require(value,message):
    if not value:raise ValueError(message)


def verify(structures,results,index):
    nodes=structures['nodes'];xyz=np.asarray(nodes['positions_m']);u=np.asarray(nodes['displacement_m']);rot=np.asarray(nodes['rotation_rad'])
    counts={}
    for family,kind,index_field,node_field in (('members','line','member_index','node_index'),
            ('plates','triangle','plate_index','structural_node_index')):
        entries=structures[family];active=[i for i,p in enumerate(entries) if p['active']]
        filename=('beams' if family=='members' else 'plates')+f'-{index:03}.vtu'
        path=results/filename
        if not active:
            # beam_io emits an empty beam frame even for a plate-only system;
            # plate_io is called only when plate definitions are present.
            if family=='plates' and not entries:
                counts[family]=0;continue
            piece=ET.parse(path).find('UnstructuredGrid/Piece')
            require(piece is not None and piece.get('NumberOfCells')=='0' and piece.get('NumberOfPoints')=='0','inactive structural VTU is not empty')
            if family=='plates':
                piece=ET.parse(results/f'plate-samples-{index:03}.vtu').find('UnstructuredGrid/Piece')
                require(piece is not None and piece.get('NumberOfCells')=='0','inactive plate samples remain visible')
            counts[family]=0;continue
        vtk=meshio.read(path);used=np.asarray(vtk.point_data[node_field],dtype=int)
        require(np.array_equal(vtk.points,xyz[used]),'structural viewer positions differ from VTU')
        require(np.array_equal(vtk.point_data['total_displacement_from_cad_m'],u[used]),'structural CAD displacement differs from VTU')
        require(np.array_equal(vtk.point_data['total_rotation_rad'],rot[used]),'structural rotation differs from VTU')
        ids=np.asarray(vtk.cell_data[index_field][0],dtype=int)
        require(set(ids)==set(active) and len(ids)==len(active),'active structural cell inventory differs from VTU')
        for row,i in enumerate(ids):
            entry=entries[i];report=entry['result']
            require(np.array_equal(used[vtk.cells_dict[kind][row]],entry['nodes']),'structural connectivity differs from VTU')
            for a,name in enumerate('xyz'):
                require(np.array_equal(vtk.cell_data['local_axis_'+name][0][row],report['local_axes'][a]),'structural local axes differ from VTU')
            if family=='members':
                require(vtk.cell_data['length_m'][0][row]==report['length_m'],'beam length differs from VTU')
                for end,start in (('i',0),('j',6)):
                    for name,begin in (('force_local_kn',start),('moment_local_kn_m',start+3)):
                        require(np.array_equal(vtk.cell_data[f'end_{end}_{name}'][0][row],report['local_end_actions'][begin:begin+3]),'beam end action differs from VTU')
                for sample,name in enumerate(('start','mid','end')):
                    for key,begin in (('force_local_kn',0),('moment_local_kn_m',3)):
                        require(np.array_equal(vtk.cell_data[f'section_{name}_{key}'][0][row],report['section_forces_at_start_mid_end'][sample][begin:begin+3]),'beam section action differs from VTU')
            else:
                for name in ('area_m2','membrane_energy_kn_m','bending_energy_kn_m','drilling_penalty_energy_kn_m'):
                    require(vtk.cell_data[name][0][row]==report[name],'plate area/energy differs from VTU')
        if family=='plates':
            raw=meshio.read(results/f'plate-samples-{index:03}.vtu');sample_keys=[]
            for row,(plate,sample) in enumerate(zip(raw.point_data['plate_index'],raw.point_data['sample_index'],strict=True)):
                plate,sample=int(plate),int(sample)
                require(plate in active and 0<=sample<3,'invalid raw plate sample identity')
                entry=entries[plate];point=entry['result']['samples'][sample];sample_keys.append((plate,sample))
                require(np.array_equal(raw.points[row],np.asarray(point['barycentric'])@xyz[entry['nodes']]),'plate sample coordinates differ from VTU')
                for name in ('membrane_strain','curvature_per_m','membrane_force_kn_m','bending_moment_kn','drilling_mismatch_rad'):
                    require(np.array_equal(raw.point_data[name][row],point[name]),'raw plate viewer sample differs from VTU: '+name)
            require(set(sample_keys)=={(p,s) for p in active for s in range(3)} and len(sample_keys)==3*len(active),'raw plate sample inventory differs')
        counts[family]=len(active)
    return counts
