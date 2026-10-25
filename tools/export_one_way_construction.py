"""Active-domain VTU, raw quadrature and separate physical/construction clocks."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import meshio
import numpy as np
import fem_reference as fem
from high_order_io import state_fields, groups, average
from one_way_reference import reconstructed_samples
from verify_acceptance_packet import sha

def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def write(path,data):
    with path.open('x',encoding='utf-8') as stream:json.dump(data,stream,indent=2)

def export(mesh_path,result_path,solver,output):
    mesh=read(mesh_path);run=read(result_path)
    assert run['complete'] and run['schema']=='soil-fem-one-way-construction-result/1'
    assert sha(solver)==run['solver_binary_sha256'];output.mkdir(parents=True,exist_ok=False)
    write(output/'identity-ground.json',run['ground']);write(output/'identity-case.json',run['case'])
    command=[str(solver),'one-way-construction-identity',str(mesh_path),str(output/'identity-ground.json'),str(output/'identity-case.json'),str(output/'identity.json')]
    log=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (output/'identity.log').write_text(log.stdout+log.stderr,encoding='utf-8');assert log.returncode==0
    assert all(read(output/'identity.json')[k]==run[k] for k in ('mesh_sha256','ground_sha256','case_sha256'))
    xyz=np.asarray(mesh['nodes']);cells=fem.cells(mesh);kind='tetra10' if cells.shape[1]==10 else 'tetra'
    layout=run['mechanical_integration_points'];reference=np.asarray(run['initialization']['initial_reference_displacement_m'])
    initial=next(f for f in run['initialization']['frames'] if f['initialization_progress']==1.)
    initial=dict(name='initial_equilibrium',phase='initial_equilibrium',progress=1.,state=initial['state'],signed_global_water_budget_error_m3=0.,cumulative_removed_water_m3=0.)
    frames=[initial]+run['frames'];rows=[];physical={}
    roots={k:ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian') for k in ('analysis-frame','physical-time','hydraulic-points','constitutive-points')}
    collections={k:ET.SubElement(root,'Collection') for k,root in roots.items()}
    for i,frame in enumerate(frames):
        s=frame['state'];active=np.asarray(s['mechanical']['active_elements'],dtype=bool);indices=np.flatnonzero(active)
        fluid=s['fluid'];hyd=fluid['hydraulic_points'];sampled=fluid['mechanical_points']
        current_ids=s['mechanical'].get('element_material_ids') or mesh['material_ids']
        if not hyd:
            current_mesh=dict(mesh,material_ids=current_ids)
            hyd,sampled=reconstructed_samples(current_mesh,dict(ground=run['ground'],case=run['case']['one_way']),s['hydraulic'])
        selected=[j for j,p in enumerate(layout) if active[p['element_index']]];points=[layout[j] for j in selected]
        assert [(p['element_index'],p['local_index']) for p in points]==[(p['element_index'],p['local_index']) for p in sampled]
        fields=('physical_pore_pressure_kpa','effective_pore_stress_kpa','bishop_chi','water_content','effective_saturation','water_equivalent_mass_per_volume')
        raw={k:np.asarray([p[k] for p in hyd]) for k in fields}
        constitutive={k:np.asarray([p[k] for p in sampled]) for k in fields}
        constitutive.update(state_fields([s['mechanical']['material_states'][j] for j in selected]))
        total=constitutive['effective_stress_compression_positive_kpa'].copy();total[:,:3]+=constitutive['effective_pore_stress_kpa'][:,None];constitutive['total_stress_compression_positive_kpa']=total
        for data,ps in ((raw,hyd),(constitutive,points)):
            data['original_element_id']=np.asarray([p['element_index'] for p in ps]);data['integration_weight_m3']=np.asarray([p['weight_m3'] for p in ps])
        def active_groups(ps):
            ids=[[] for _ in cells]
            for j,p in enumerate(ps):ids[p['element_index']].append(j)
            assert all(ids[e] for e in indices)
            return ids
        hids=active_groups(hyd);mids=active_groups(points)
        # Only active cells are averaged; no empty/inactive GP groups are used.
        def means(values,ps,ids):
            return np.asarray([np.average(values[ids[e]],axis=0,weights=np.asarray([ps[j]['weight_m3'] for j in ids[e]])) for e in indices])
        cell={'original_element_id':[indices],'material_id':[np.asarray(current_ids)[active]],
              'cell_total_head_m':[run['case']['one_way']['flow'].get('head_datum_m',0.)+np.asarray([s['hydraulic']['relative_cell_head_m'][e] for e in indices])]}
        for k in fields:cell['hydraulic_volume_mean_'+k]=[means(raw[k],hyd,hids)]
        for k in ('effective_stress_compression_positive_kpa','total_stress_compression_positive_kpa','void_ratio'):
            cell['constitutive_volume_mean_'+k]=[means(constitutive[k],points,mids)]
        live=np.zeros(len(xyz),dtype=np.uint8);live[np.unique(cells[active])]=1
        u=np.asarray(s['mechanical']['displacement_m']);point={'original_node_id':np.arange(len(xyz)),'node_active':live,'displacement_from_seed_m':u,'displacement_from_initial_equilibrium_m':u-reference}
        name=f'soil-{i:04}.vtu';meshio.write_points_cells(output/name,xyz,[(kind,cells[active])],point_data=point,cell_data=cell)
        check=meshio.read(output/name);assert np.array_equal(check.cells_dict[kind],cells[active])
        for k,v in point.items():assert np.array_equal(check.point_data[k],v)
        for k,v in cell.items():assert np.array_equal(check.cell_data[k][0],v[0])
        ET.SubElement(collections['analysis-frame'],'DataSet',timestep=str(i),group='',part='0',file=name)
        physical[s['hydraulic']['time_s']]=name
        for label,ps,data in (('hydraulic-points',hyd,raw),('constitutive-points',points,constitutive)):
            name=f'{label}-{i:04}.vtu';positions=np.asarray([p['position_m'] for p in ps]);meshio.write_points_cells(output/name,positions,[('vertex',np.arange(len(ps))[:,None])],point_data=data)
            check=meshio.read(output/name);assert np.array_equal(check.points,positions)
            for k,v in data.items():assert np.array_equal(check.point_data[k],v,equal_nan=True)
            ET.SubElement(collections[label],'DataSet',timestep=str(i),group='',part='0',file=name)
        rows.append(dict(analysis_frame=i,name=frame['name'],phase=frame['phase'],physical_time_s=s['hydraulic']['time_s'],release_progress=frame['progress'],active_elements=len(indices),water_m3=fluid['water_equivalent_mass_volume_m3'],removed_water_m3=frame['cumulative_removed_water_m3'],added_water_m3=frame.get('cumulative_added_water_m3',0.),water_budget_error_m3=frame['signed_global_water_budget_error_m3']))
    for time,name in sorted(physical.items()):ET.SubElement(collections['physical-time'],'DataSet',timestep=str(time),group='',part='0',file=name)
    for key,root in roots.items():ET.ElementTree(root).write(output/f'{key}.pvd',encoding='utf-8',xml_declaration=True)
    for key in roots:
        entries=ET.parse(output/f'{key}.pvd').find('Collection');times=[float(row.attrib['timestep']) for row in entries]
        assert all(b>a for a,b in zip(times,times[1:]));assert all((output/row.attrib['file']).exists() for row in entries)
    for i,record in enumerate(run['excavations']):
        archived=record.get('removed_material_points',[])
        ids=record['removed_mechanical_points'];data=state_fields([p['state'] for p in archived] if archived else [run['last_committed']['mechanical']['material_states'][j] for j in ids])
        if archived:
            assert ids==[p['integration_point_index'] for p in archived]
            data['archived_material_id']=np.asarray([p['material_id'] for p in archived]);data['active_before']=np.asarray([p['active_before'] for p in archived],dtype=np.uint8)
        data['original_integration_point_id']=np.asarray(ids);positions=np.asarray([layout[j]['position_m'] for j in ids])
        path=output/f'removed-history-{i:03}.vtu';meshio.write_points_cells(path,positions,[('vertex',np.arange(len(ids))[:,None])],point_data=data)
        check=meshio.read(path)
        for k,v in data.items():assert np.array_equal(check.point_data[k],v,equal_nan=True)
    for i,record in enumerate(run.get('transitions',[])):
        ps=record['old_material_points'];ids=[p['integration_point_index'] for p in ps];positions=np.asarray([layout[j]['position_m'] for j in ids])
        data=state_fields([p['state'] for p in ps]);data['original_integration_point_id']=np.asarray(ids);data['archived_material_id']=np.asarray([p['material_id'] for p in ps]);data['active_before']=np.asarray([p['active_before'] for p in ps],dtype=np.uint8)
        path=output/f'transition-old-history-{i:03}.vtu';meshio.write_points_cells(path,positions,[('vertex',np.arange(len(ids))[:,None])],point_data=data);check=meshio.read(path)
        for k,v in data.items():assert np.array_equal(check.point_data[k],v,equal_nan=True)
    with (output/'analysis-frames.csv').open('x',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    write(output/'export-verification.json',dict(schema='soil-fem-one-way-construction-export/1',passed=True,analysis_frames=len(frames),physical_frames=len(physical),all_vtu_and_pvd_roundtrip='passed',transfer_point_source='stored' if run['case']['one_way']['store_transfer_integration_points'] else 'reconstructed',mesh_sha256=sha(mesh_path),result_sha256=sha(result_path),solver_sha256=sha(solver),scope='Analysis-frame PVD uses frame indices, not seconds. Physical-time PVD selects the last accepted state at each time, including final excavation release. Raw active GP layouts and removed histories are separate; no nodal water-pressure smoothing.'))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mesh',type=Path);parser.add_argument('result',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--solver',type=Path,required=True);a=parser.parse_args();export(a.mesh.resolve(),a.result.resolve(),a.solver.resolve(),a.output.resolve())
