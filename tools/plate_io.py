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
    offset_scope=any(p.get('formulation')=='mindlin_tri6_offset' for p in cfg['plates'])
    roots=[ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian') for _ in range(2)]
    collections=[ET.SubElement(r,'Collection') for r in roots];rows=[]
    for fi,frame in enumerate(run['frames']):
        report=frame['beam_results'];state=frame['state']['beams'];active=report.get('plates',[])
        assert len(active)==sum(state['active_plates'])
        filename=f'plates-{fi:03}.vtu';pointfile=f'plate-samples-{fi:03}.vtu'
        if not active:
            empty_frame(output/filename);empty_frame(output/pointfile)
        else:
            ids=np.array([plates[p['id']] for p in active]);cells=[[nodes[n] for n in cfg['plates'][i]['nodes']] for i in ids]
            used=np.unique([n for cell in cells for n in cell]);mapping=np.full(len(nodes),-1);mapping[used]=np.arange(len(used))
            xyz=np.asarray(report['node_positions_m']);u=np.asarray(state['node_displacement_m']);rotation=np.asarray(state['node_rotation_rad'])
            points=dict(structural_node_index=used,total_displacement_from_cad_m=u[used],total_rotation_rad=rotation[used])
            fields=dict(plate_index=ids,area_m2=np.array([p['area_m2'] for p in active]))
            if offset_scope:
                fields['attachment_offset_m']=np.array([cfg['plates'][i].get('attachment_offset_m',0.) for i in ids])
                neutral_xyz=xyz.copy();neutral_u=u.copy()
                for i,p,cell in zip(ids,active,cells,strict=True):
                    arm=np.asarray(p['local_axes'][2])*cfg['plates'][i].get('attachment_offset_m',0.)
                    neutral_xyz[cell]=xyz[cell]-arm
                    neutral_u[cell]=u[cell]-np.cross(rotation[cell],arm)
                points.update(neutral_reference_position_m=neutral_xyz[used],neutral_total_displacement_m=neutral_u[used])
            for key in ('membrane_energy_kn_m','bending_energy_kn_m','drilling_penalty_energy_kn_m'):
                fields[key]=np.array([p[key] for p in active])
            for a,name in enumerate('xyz'):fields[f'local_axis_{name}']=np.array([p['local_axes'][a] for p in active])
            # Preserve each kind as its native VTK cell, including mixed inventories.
            groups=[(kind,np.array([i for i,c in enumerate(cells) if len(c)==count],dtype=int)) for kind,count in (('triangle',3),('triangle6',6))]
            groups=[(kind,indices) for kind,indices in groups if len(indices)]
            blocks=[(kind,mapping[np.array([cells[i] for i in indices])]) for kind,indices in groups]
            if any(len(c)==6 for c in cells):fields['transverse_shear_energy_kn_m']=np.array([p.get('transverse_shear_energy_kn_m',0.) for p in active])
            meshio.write_points_cells(output/filename,xyz[used],blocks,point_data=points,cell_data={k:[v[indices] for _,indices in groups] for k,v in fields.items()})
            check=meshio.read(output/filename)
            for kind,block in blocks:assert np.array_equal(check.cells_dict[kind],block)
            for k,v in points.items():assert np.array_equal(check.point_data[k],v)
            for k,v in fields.items():
                for block,(_,indices) in zip(check.cell_data[k],groups,strict=True):assert np.array_equal(block,v[indices])
            positions=[];raw={k:[] for k in ('plate_index','sample_index','membrane_strain','curvature_per_m','membrane_force_kn_m','bending_moment_kn','drilling_mismatch_rad')}
            high=any(len(c)==6 for c in cells)
            if high:raw.update(transverse_shear_strain=[],transverse_shear_force_kn_m=[])
            for plate_id,p,cell in zip(ids,active,cells,strict=True):
                assert len(p['samples'])==(6 if len(cell)==6 else 3)
                for qi,sample in enumerate(p['samples']):
                    bary=np.asarray(sample['barycentric']);assert np.all(bary>=0.) and abs(bary.sum()-1.)<1e-14
                    arm=np.asarray(p['local_axes'][2])*cfg['plates'][plate_id].get('attachment_offset_m',0.) if offset_scope else np.zeros(3)
                    positions.append(bary@xyz[cell[:3]]-arm);raw['plate_index'].append(plate_id);raw['sample_index'].append(qi)
                    row=dict(stage_index=fi,stage=frame['name'],plate_id=p['id'],sample_index=qi)
                    for key in raw:
                        if key not in ('plate_index','sample_index'):raw[key].append(sample.get(key,[0.,0.]) if key.startswith('transverse_shear') else sample[key])
                    for key in ('barycentric','membrane_strain','curvature_per_m','membrane_force_kn_m','bending_moment_kn'):
                        for a,value in enumerate(sample[key]):row[f'{key}_{a}']=value
                    if high:
                        for key in ('transverse_shear_strain','transverse_shear_force_kn_m'):
                            for a,value in enumerate(sample.get(key,[0.,0.])):row[f'{key}_{a}']=value
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
        writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for row in rows for k in row)) if rows else ['stage_index','plate_id','sample_index'])
        writer.writeheader();writer.writerows(rows)
    metadata=dict(schema='soil-fem-plate-export/1',plates=plates,nodes=nodes,
        geometry='Undeformed flat CAD triangles; node displacement is from CAD, separate from soil display RESET.',
        fields='Raw three-point DKT samples in local axes. Membrane force kN/m; bending moment kN m/m; curvature 1/m. No smoothing.',
        drilling='Explicit membrane-spin penalty, with energy reported separately from membrane and bending energy.',
        roundtrip='Nonempty triangle and vertex VTU data read back exactly; empty active sets checked as zero-cell VTK XML.')
    if any(p.get('formulation')=='mindlin_tri6' for p in cfg['plates']):
        metadata.update(schema='soil-fem-plate-export/2',fields='Raw DKT (3 points) or Mindlin Tri6 (6 points) samples in local axes; transverse shear strain and force kN/m; separate shear energy. No smoothing.',geometry='Undeformed same-reference-plane triangle/triangle6 cells; offsets unsupported.')
    if offset_scope:
        metadata.update(schema='soil-fem-plate-export/3',geometry='triangle6 points and total displacement refer to the attachment surface; neutral reference position and displacement are explicit point fields; raw sample positions lie on the neutral plane.',fields='Six raw Mindlin points in local axes, neutral-plane strain/curvature/resultants and four energies; attachment_offset_m is signed neutral-to-surface distance.',native_gui_accepted=False)
    (output/'plate-fields.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
