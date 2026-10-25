"""Raw JOINT quadrature histories; inactive points remain explicitly archived."""
import csv
import xml.etree.ElementTree as ET
import meshio
import numpy as np


def export(result,output):
    layout=result['joint_integration_points'];xyz=np.asarray([p['position_m'] for p in layout]);weights=np.asarray([p['weight_m2'] for p in layout])
    root=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian');collection=ET.SubElement(root,'Collection');rows=[]
    for fi,frame in enumerate(result['frames']):
        state=frame['state'];points=state['joints']['points'];report=frame['joint_results'];u=np.asarray(state['displacement_m'])
        reference=np.asarray(state.get('display_reference_displacement_m') or result['reference_displacement_m']);fields={}
        for key in ('normal','jump_m','traction_kpa','plastic_slip_m','accumulated_plastic_slip_m','plastic_dissipation_kpa_m','numerical_dissipation_kpa_m','endpoint_work_kpa_m'):fields[key]=np.asarray([s[key] for s in points])
        for key in ('active','gap_m','opening_m','normal_pressure_kpa','shear_traction_kpa','elastic_energy_kpa_m','contact_reference_offset_m'):fields[key]=np.asarray([p[key] for p in report['points']],dtype=int if key=='active' else float)
        if report.get('water') is not None:
            for key in ('physical_pore_pressure_kpa','effective_contact_pressure_kpa','total_normal_pressure_kpa','friction_limit_kpa'):fields[key]=np.asarray([p['hydraulic'][key] for p in report['points']])
            fields['pressure_source']=np.asarray([{'soil':0,'continuous':1,'prescribed':2}[p['hydraulic']['source']] for p in report['points']])
        fields['status']=np.asarray([{'open':0,'stick':1,'slip':2}[s['status']] for s in points])
        fields['integration_point_id']=np.arange(len(points));fields['integration_weight_m2']=weights
        fields['owner_a']=np.asarray([p['owner_a'] for p in layout]);fields['owner_b']=np.asarray([p['owner_b'] for p in layout])
        for side in ('a','b'):fields['side_'+side+'_total_displacement_m']=np.asarray([np.asarray(p['shape'])@u[p['nodes_'+side]] for p in layout])
        fields['displacement_m']=np.asarray([np.asarray(p['shape'])@((u-reference)[p['nodes_a']]+(u-reference)[p['nodes_b']])/2 for p in layout])
        name=f'joints-{fi:03}.vtu';meshio.write_points_cells(output/name,xyz,[('vertex',np.arange(len(points))[:,None])],point_data=fields)
        checked=meshio.read(output/name)
        assert np.array_equal(checked.points,xyz)
        for key,value in fields.items():assert np.array_equal(checked.point_data[key],value),key
        ET.SubElement(collection,'DataSet',timestep=str(fi),group=frame['name'],part='0',file=name)
        row={'stage':frame['name'],'active_points':int(fields['active'].sum()),'open_points':int(((fields['status']==0)&(fields['active']==1)).sum())}
        row.update({k:report[k] for k in ('elastic_energy_kn_m','plastic_dissipation_kn_m','numerical_dissipation_kn_m','endpoint_work_kn_m')});rows.append(row)
        if report.get('water') is not None:row['maximum_physical_pore_pressure_kpa']=float(fields['physical_pore_pressure_kpa'].max())
    ET.ElementTree(root).write(output/'joints.pvd',encoding='utf-8',xml_declaration=True)
    with (output/'joints.csv').open('x',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    if any(f['joint_results'].get('water') is not None for f in result['frames']):
        import json
        with (output/'joint-water-fields.json').open('x',encoding='utf-8') as stream:json.dump(dict(pressure_source_codes={'0':'soil','1':'continuous','2':'prescribed'},material_stress_reference_codes={'0':'soil effective stress','1':'impermeable structure total/material stress'},scope='Effective contact spring stress and physical pressure are separate. Total normal contact pressure equals their sum; friction uses effective contact pressure. Water force integration is stored separately from GP display samples.'),stream,indent=2)
