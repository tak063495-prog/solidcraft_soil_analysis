"""Separate flat-plate cells and raw DKT samples; no nodal stress smoothing."""
import csv
import json
import xml.etree.ElementTree as ET
import meshio
import numpy as np
from beam_io import empty_frame


def export(run,output):
    cfg=run['case']['beam_system'];nodes={n['id']:i for i,n in enumerate(cfg['nodes'])}
    plates={p['id']:i for i,p in enumerate(cfg['plates'])}
    roots=[ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian') for _ in range(2)]
    collections=[ET.SubElement(r,'Collection') for r in roots];rows=[]
    for fi,frame in enumerate(run['frames']):
        report=frame['beam_results'];state=frame['state']['beams'];active=report.get('plates',[])
        assert len(active)==sum(state['active_plates'])
        filename=f'plates-{fi:03}.vtu';pointfile=f'plate-samples-{fi:03}.vtu'
        if not active:
            empty_frame(output/filename);empty_frame(output/pointfile)
        else:
            ids=np.array([plates[p['id']] for p in active]);cells=np.array([[nodes[n] for n in cfg['plates'][i]['nodes']] for i in ids])
            used=np.unique(cells);mapping=np.full(len(nodes),-1);mapping[used]=np.arange(len(used))
            xyz=np.asarray(report['node_positions_m']);u=np.asarray(state['node_displacement_m']);rotation=np.asarray(state['node_rotation_rad'])
            points=dict(structural_node_index=used,total_displacement_from_cad_m=u[used],total_rotation_rad=rotation[used])
            fields=dict(plate_index=ids,area_m2=np.array([p['area_m2'] for p in active]))
            for key in ('membrane_energy_kn_m','bending_energy_kn_m','drilling_penalty_energy_kn_m'):
                fields[key]=np.array([p[key] for p in active])
            for a,name in enumerate('xyz'):fields[f'local_axis_{name}']=np.array([p['local_axes'][a] for p in active])
            meshio.write_points_cells(output/filename,xyz[used],[('triangle',mapping[cells])],point_data=points,cell_data={k:[v] for k,v in fields.items()})
            check=meshio.read(output/filename);assert np.array_equal(check.cells_dict['triangle'],mapping[cells])
            for k,v in points.items():assert np.array_equal(check.point_data[k],v)
            for k,v in fields.items():assert np.array_equal(check.cell_data[k][0],v)
            positions=[];raw={k:[] for k in ('plate_index','sample_index','membrane_strain','curvature_per_m','membrane_force_kn_m','bending_moment_kn','drilling_mismatch_rad')}
            for plate_id,p,cell in zip(ids,active,cells,strict=True):
                assert len(p['samples'])==3
                for qi,sample in enumerate(p['samples']):
                    bary=np.asarray(sample['barycentric']);assert np.all(bary>=0.) and abs(bary.sum()-1.)<1e-14
                    positions.append(bary@xyz[cell]);raw['plate_index'].append(plate_id);raw['sample_index'].append(qi)
                    row=dict(stage_index=fi,stage=frame['name'],plate_id=p['id'],sample_index=qi)
                    for key in raw:
                        if key not in ('plate_index','sample_index'):raw[key].append(sample[key])
                    for key in ('barycentric','membrane_strain','curvature_per_m','membrane_force_kn_m','bending_moment_kn'):
                        for a,value in enumerate(sample[key]):row[f'{key}_{a}']=value
                    row['drilling_mismatch_rad']=sample['drilling_mismatch_rad'];rows.append(row)
            raw={k:np.asarray(v) for k,v in raw.items()}
            meshio.write_points_cells(output/pointfile,np.asarray(positions),[('vertex',np.arange(len(positions))[:,None])],point_data=raw)
            checked=meshio.read(output/pointfile)
            for k,v in raw.items():assert np.array_equal(checked.point_data[k],v)
        for collection,name in zip(collections,(filename,pointfile),strict=True):
            ET.SubElement(collection,'DataSet',timestep=str(fi),group=frame['name'],part='0',file=name)
    for root,name in zip(roots,('plates.pvd','plate-samples.pvd'),strict=True):
        ET.ElementTree(root).write(output/name,encoding='utf-8',xml_declaration=True)
    with (output/'plate-samples.csv').open('x',encoding='utf-8-sig',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]) if rows else ['stage_index','plate_id','sample_index'])
        writer.writeheader();writer.writerows(rows)
    metadata=dict(schema='soil-fem-plate-export/1',plates=plates,nodes=nodes,
        geometry='Undeformed flat CAD triangles; node displacement is from CAD, separate from soil display RESET.',
        fields='Raw three-point DKT samples in local axes. Membrane force kN/m; bending moment kN m/m; curvature 1/m. No smoothing.',
        drilling='Explicit membrane-spin penalty, with energy reported separately from membrane and bending energy.',
        roundtrip='Nonempty triangle and vertex VTU data read back exactly; empty active sets checked as zero-cell VTK XML.')
    (output/'plate-fields.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
