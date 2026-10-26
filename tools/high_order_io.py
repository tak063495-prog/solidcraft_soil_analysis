"""VTU/PVD: keep volume averages distinct from raw integration-point fields."""
import hashlib
import json
from pathlib import Path
import csv
import xml.etree.ElementTree as ET
import meshio
import numpy as np


def point_state(point):
    return point.get('mcc_state',point.get('mohr_coulomb_state',point.get('invariant_double_hardening_state')))


def state_fields(states):
    fields={name:np.asarray([s[name] for s in states]) for name in
        ('effective_stress_compression_positive_kpa','plastic_strain_compression_positive','void_ratio')}
    fields['material_model']=np.asarray([0 if 'preconsolidation_pressure_kpa' in s else (2 if 'cap_pressure_kpa' in s else 1) for s in states])
    for name in ('preconsolidation_pressure_kpa','equivalent_plastic_strain','plastic_work_kpa',
                 'cap_pressure_kpa','initial_cap_pressure_kpa','shear_plastic_strain','cap_plastic_volumetric_strain'):
        if any(name in s for s in states):fields[name]=np.asarray([s.get(name,np.nan) for s in states])
    return fields


def groups(layout,count):
    ids=[[] for _ in range(count)]
    for i,p in enumerate(layout):ids[p['element_index']].append(i)
    assert all(ids)
    return ids


def average(values,layout,ids):
    v=np.asarray(values);weights=np.asarray([p['weight_m3'] for p in layout])
    return np.asarray([np.average(v[g],axis=0,weights=weights[g]) for g in ids])


def initial(mesh,state,output):
    output=Path(output);point_output=output.with_name(output.stem+'-points.vtu')
    if output.exists() or point_output.exists():raise ValueError('initial VTU output exists')
    points=state['point_states'];layout=[{'element_index':p['element_index'],'weight_m3':p['weight_m3']} for p in points]
    ids=groups(layout,len(mesh['tetrahedra']));fields={}
    for name in ('pore_pressure_kpa','porosity','mcc_pc_over_mean_effective','normalized_yield','bulk_modulus_kpa','shear_modulus_kpa','total_stress_compression_positive_kpa'):
        fields[name]=np.asarray([p[name] for p in points],dtype=float)
    fields.update(state_fields([point_state(p) for p in points]))
    cell={'material_id':[np.asarray(mesh['material_ids'])],'region_id':[np.asarray(mesh['region_ids'])],'integration_point_count':[np.array([len(g) for g in ids])]}
    cell.update({'volume_average_'+name:[average(v,layout,ids)] for name,v in fields.items()})
    quadratic=bool(mesh.get('tetrahedra10'));kind='tetra10' if quadratic else 'tetra';cells=mesh['tetrahedra10'] if quadratic else mesh['tetrahedra']
    meshio.write_points_cells(output,mesh['nodes'],[(kind,cells)],cell_data=cell)
    meshio.write_points_cells(point_output,[p['centroid_m'] for p in points],[('vertex',np.arange(len(points))[:,None])],
        point_data={**fields,'original_element_id':np.array([p['element_index'] for p in points]),'integration_weight_m3':np.array([p['weight_m3'] for p in points])})
    checked=meshio.read(output);checked_points=meshio.read(point_output)
    assert np.array_equal(checked.cells_dict[kind],cells)
    assert np.allclose(checked_points.point_data['void_ratio'],fields['void_ratio'])


def stages(mesh_file,result_file,output):
    fingerprints={str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (mesh_file,result_file)}
    mesh=json.loads(Path(mesh_file).read_text(encoding='utf-8-sig'));result=json.loads(Path(result_file).read_text(encoding='utf-8-sig'))
    mixed=result['schema'] in ('soil-fem-mixed-mcc-stage-result/1','soil-fem-projected-mcc-stage-result/1','soil-fem-mixed-mcc-stage-result/2','soil-fem-projected-mcc-stage-result/2')
    if not mixed and result['schema'] not in ('soil-fem-stage-result/2','soil-fem-stage-result/3','soil-fem-stage-result/4','soil-fem-stage-result/5','soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/13','soil-fem-stage-result/14') or not result['complete']:raise ValueError('completed high-order/general-material result required')
    if mixed:
        import mixed_mcc_io
        mixed_mcc_io.validate(mesh,result)
    elif any(f['state'].get('mixed_mcc') is not None for f in result['frames']):
        raise ValueError('mixed primary state requires its own result schema')
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    quadratic=bool(mesh.get('tetrahedra10'));kind='tetra10' if quadratic else 'tetra'
    xyz=np.array(mesh['nodes']);cells=np.array(mesh['tetrahedra10'] if quadratic else mesh['tetrahedra']);layout=result['integration_points'];ids=groups(layout,len(cells))
    initial_reference=np.asarray(result.get('reference_displacement_m',np.zeros_like(xyz)))
    roots=[ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian') for _ in range(2)]
    collections=[ET.SubElement(r,'Collection') for r in roots];rows=[]
    for frame_index,frame in enumerate(result['frames']):
        state=frame['state'];reference=np.asarray(state.get('display_reference_displacement_m') or initial_reference);active=np.asarray(state['active_elements']);used=np.unique(cells[active]);mapping=np.full(len(xyz),-1);mapping[used]=np.arange(len(used))
        u=np.asarray(state['displacement_m']);old=np.asarray(result['frames'][max(0,frame_index-1)]['state']['displacement_m'])
        if frame_index and result['case']['stages'][frame_index-1]['kind']=='birth':old=np.asarray(next(a for a in result['transition_audits'] if a['stage']==frame['name'])['prepared_state']['displacement_m'])
        raw=state_fields(state['material_states'])
        mixed_fields=None
        if mixed:
            mixed_fields=mixed_mcc_io.frame_fields(mesh,result,frame)
            raw.update(mixed_fields['raw'])
            mixed_mcc_io.write_primary(output,frame_index,mesh,mixed_fields)
            mixed_mcc_io.write_history(output,frame_index,frame,mixed_fields)
        h=state.get('hydrostatic_water_level_m',result['ground']['initialization']['water_level_m']);gamma=result['ground']['initialization']['water_unit_weight_kn_m3']
        pore=np.zeros(len(layout)) if h is None else np.maximum(h-np.array([p['position_m'][2] for p in layout]),0)*gamma
        if state.get('pore_pressure_kpa') is not None:pore=np.asarray(state['pore_pressure_kpa'])
        structures=set((result['case'].get('joint_system') or {}).get('structure_bodies',[]))
        if result['schema']=='soil-fem-stage-result/13':
            bodies={r['id']:r['body'] for r in mesh['regions']};total=np.asarray([bodies[mesh['region_ids'][p['element_index']]] in structures for p in layout])
            if state.get('pore_pressure_kpa') is not None:assert np.all(pore[total]==0.)
            pore=pore.copy();pore[total]=0.;raw['material_stress_reference']=total.astype(np.uint8)
        raw['pore_pressure_kpa']=pore
        raw['total_stress_compression_positive_kpa']=raw['effective_stress_compression_positive_kpa'].copy()
        raw['total_stress_compression_positive_kpa'][:,:3]+=pore[:,None]
        assert all(len(v)==len(layout) for v in raw.values())
        cell={'original_element_id':[np.flatnonzero(active)],'material_id':[np.asarray(state.get('element_material_ids') or mesh['material_ids'])[active]],'region_id':[np.asarray(mesh['region_ids'])[active]],
              'integration_point_count':[np.array([len(g) for g in ids])[active]]}
        cell.update({'volume_average_'+name:[average(value,layout,ids)[active]] for name,value in raw.items()})
        filename=f'stage-{frame_index:03}.vtu'
        node_fields={'original_node_id':used,'displacement_m':(u-reference)[used],'stage_displacement_increment_m':(u-old)[used]}
        if mixed or result['schema'] in ('soil-fem-stage-result/4','soil-fem-stage-result/5','soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/13','soil-fem-stage-result/14'):node_fields.update(total_displacement_from_seed_m=u[used],initialization_displacement_m=initial_reference[used])
        if mixed or result['schema'] in ('soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/13','soil-fem-stage-result/14'):node_fields['display_reference_displacement_m']=reference[used]
        meshio.write_points_cells(output/filename,xyz[used],[(kind,mapping[cells[active]])],point_data=node_fields,cell_data=cell)
        qp=np.array([i for i,p in enumerate(layout) if active[p['element_index']]])
        positions=np.array([layout[i]['position_m'] for i in qp])
        point_u=np.array([np.array(layout[i]['shape_functions'])@(u-reference)[cells[layout[i]['element_index']]] for i in qp])
        point_file=f'integration-points-{frame_index:03}.vtu'
        point_fields={**{k:v[qp] for k,v in raw.items()},'integration_point_id':qp,
                'original_element_id':np.array([layout[i]['element_index'] for i in qp]),
                'integration_weight_m3':np.array([layout[i]['weight_m3'] for i in qp]),'displacement_m':point_u}
        if mixed or result['schema'] in ('soil-fem-stage-result/5','soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/13','soil-fem-stage-result/14'):point_fields['material_id']=np.array([(state.get('element_material_ids') or mesh['material_ids'])[layout[i]['element_index']] for i in qp])
        if mixed or result['schema'] in ('soil-fem-stage-result/4','soil-fem-stage-result/5','soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/13','soil-fem-stage-result/14'):
            point_fields['total_displacement_from_seed_m']=np.array([np.array(layout[i]['shape_functions'])@u[cells[layout[i]['element_index']]] for i in qp])
            point_fields['initialization_displacement_m']=np.array([np.array(layout[i]['shape_functions'])@initial_reference[cells[layout[i]['element_index']]] for i in qp])
        meshio.write_points_cells(output/point_file,positions,[('vertex',np.arange(len(qp))[:,None])],point_data=point_fields)
        checked=meshio.read(output/filename);point_checked=meshio.read(output/point_file)
        assert np.array_equal(checked.cells_dict[kind],mapping[cells[active]])
        assert np.array_equal(checked.cell_data['material_id'][0],np.asarray(state.get('element_material_ids') or mesh['material_ids'])[active])
        if mixed or result['schema'] in ('soil-fem-stage-result/5','soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/13','soil-fem-stage-result/14'):assert np.array_equal(point_checked.point_data['material_id'],point_fields['material_id'])
        assert np.allclose(checked.point_data['displacement_m'],(u-reference)[used])
        for k,v in raw.items():assert np.allclose(point_checked.point_data[k],v[qp],equal_nan=True)
        if mixed:
            for k,v in raw.items():assert np.array_equal(point_checked.point_data[k],v[qp],equal_nan=True)
            for k,v in cell.items():assert np.array_equal(checked.cell_data[k][0],v[0],equal_nan=True)
        assert np.allclose(point_checked.point_data['displacement_m'],point_u)
        for collection,file in zip(collections,(filename,point_file)):
            ET.SubElement(collection,'DataSet',timestep=str(frame_index),group=frame['name'],part='0',file=file)
        row={'stage_index':frame_index,'stage':frame['name'],'active_elements':int(active.sum()),'active_material_points':len(qp),'relative_force_residual':frame['relative_force_residual']}
        if mixed:row.update(mixed_fields['metrics'])
        if 'hydrostatic_water_level_m' in state:row['hydrostatic_water_level_m']=h
        for name,m in frame['monitors'].items():
            row[name+'_u_z_m']=m['mean_displacement_z_m'];row[name+'_increment_z_m']=m['stage_increment_z_m']
            row[name+'_minimum_u_z_m']=m['minimum_displacement_z_m'];row[name+'_maximum_u_z_m']=m['maximum_displacement_z_m']
        rows.append(row)
    for root,name in zip(roots,('stages.pvd','integration-points.pvd')):ET.ElementTree(root).write(output/name,encoding='utf-8',xml_declaration=True)
    for name,data in [('stages.csv',rows),('attempts.csv',result['attempts'])]:
        with (output/name).open('x',encoding='utf-8-sig',newline='') as f:
            fields=list(dict.fromkeys(k for row in data for k in row)) if data else ['stage','accepted'];writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(data)
    if result.get('geostatic') is not None:
        data=result['geostatic']['attempts']
        with (output/'geostatic-attempts.csv').open('x',encoding='utf-8-sig',newline='') as f:
            fields=list(dict.fromkeys(k for row in data for k in row)) if data else ['stage','accepted'];writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(data)
        g=result['geostatic'];summary={k:g[k] for k in ('complete','progress','maximum_absolute_engineering_strain_from_seed','final_relative_physical_force_residual')}
        summary['accepted_increments']=sum(a['accepted'] for a in data);summary['rejected_increments']=sum(not a['accepted'] for a in data)
        (output/'geostatic-summary.json').write_text(json.dumps(summary,indent=2))
    if result['case'].get('joint_system'):
        from joint_io import export
        export(result,output)
    if result['case'].get('beam_system'):
        from beam_io import export
        export(result,output)
        if result['case']['beam_system'].get('plates'):
            from plate_io import export as export_plates
            export_plates(result,output)
    if result['case'].get('anchor_system'):
        from anchor_io import export
        export(result,output)
    if result.get('flow_audits'):
        from export_darcy import export
        for audit in result['flow_audits']:
            hydraulic_mesh=audit['joint_flow']['hydraulic_mesh'] if audit.get('joint_flow') else mesh
            export(hydraulic_mesh,audit['solution'],output/f"flow-{audit['stage_index']:03}")
            if audit.get('joint_flow'):
                with (output/f"joint-flow-map-{audit['stage_index']:03}.json").open('x',encoding='utf-8') as stream:json.dump({k:audit['joint_flow'][k] for k in ('original_mesh_sha256','original_case_sha256','original_node_to_hydraulic_node','original_field')},stream,indent=2)
    assert fingerprints=={str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (mesh_file,result_file)},'export inputs changed'
    record={'schema':'soil-fem-stage-export/2','result_sha256':hashlib.sha256(Path(result_file).read_bytes()).hexdigest(),
            'mesh_sha256':hashlib.sha256(Path(mesh_file).read_bytes()).hexdigest(),'vtu_roundtrip':'passed',
            'integration_point_roundtrip':'passed','cell_field_meaning':'volume_average','frame_count':len(rows)}
    if result['case'].get('joint_system'):record['joint_point_roundtrip']='passed'
    if (result['case'].get('beam_system') or {}).get('plates'):record['plate_sample_roundtrip']='passed'
    if mixed:
        mixed_mcc_io.write_collection(output,result['frames'])
        record.update(schema='soil-fem-stage-export/3',formulation='Tet10 / material-split P1 effective mean pressure / GP volume MCC',
            tangent_method=result['frames'][0]['state']['mixed_mcc'].get('tangent_method',mixed_mcc_io.LEGACY_TANGENT),
            pressure_primary_roundtrip='passed',all_gp_history_csv_roundtrip='passed',
            shared_corner_material_pressures_averaged=False,material_and_kinematic_volume_kept_distinct=True,
            independent_field_accuracy_checked_by_exporter=False,
            pressure_field_meaning='compression-positive mean effective stress, distinct from pore pressure',
            inactive_primary_pressure_meaning='undefined; CSV blank, frozen material history retained',
            source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(mixed_mcc_io.__file__),Path(mixed_mcc_io.ref.__file__))},
            sha256={str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in output.rglob('*') if p.is_file()})
        if result['schema'] in mixed_mcc_io.PROJECTED_SCHEMAS:
            record.update(schema='soil-fem-stage-export/4',formulation=mixed_mcc_io.PROJECTED_FORMULATION,
                pressure_consistency='weak material-split P1 projection; material GP pressure is a distinct field',
                equilibrium_stress_field='work_conjugate_effective_stress_compression_positive_kpa',
                domain_change_history='frozen integral compatibility offsets; material pc/e/plastic histories retained')
    (output/'export-verification.json').write_text(json.dumps(record,indent=2))
    print(f'Exported {len(rows)} {kind} frames and raw integration-point frames: {output}')
