"""Displacement, physical u, chi*u and water mass, with distinct GP layouts."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import meshio
import numpy as np
import fem_reference as fem
from high_order_io import state_fields,groups,average
from verify_acceptance_packet import sha
from one_way_reference import reconstructed_samples


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def export(mesh_path,result_path,solver,output):
    mesh=read(mesh_path);run=read(result_path);assert run['complete'] and run['schema']=='soil-fem-one-way-result/1'
    assert sha(solver)==run['solver_binary_sha256'];output.mkdir(parents=True,exist_ok=False)
    for name,data in [('ground',run['ground']),('case',run['case'])]:
        with (output/f'identity-{name}.json').open('x',encoding='utf-8') as f:json.dump(data,f)
    command=[str(solver),'one-way-identity',str(mesh_path),str(output/'identity-ground.json'),str(output/'identity-case.json'),str(output/'identity.json')]
    log=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace');(output/'identity.log').write_text(log.stdout+log.stderr,encoding='utf-8');assert log.returncode==0
    identity=read(output/'identity.json');assert all(identity[k]==run[k] for k in ('mesh_sha256','ground_sha256','case_sha256'))
    xyz=np.asarray(mesh['nodes']);cells=fem.cells(mesh);kind='tetra10' if cells.shape[1]==10 else 'tetra';mechanical=run['mechanical_integration_points'];mids=groups(mechanical,len(cells));rows=[]
    roots={k:ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian') for k in ('soil','hydraulic-points','constitutive-points')};collections={k:ET.SubElement(v,'Collection') for k,v in roots.items()}
    initial=next(f for f in run['frames'] if f['initialization_progress']==1.);frames=[initial]+[f for f in run['frames'] if f['phase']=='transient']
    assert all(frames[i+1]['state']['hydraulic']['time_s']>frames[i]['state']['hydraulic']['time_s'] for i in range(len(frames)-1))
    reference=np.asarray(run['initial_reference_displacement_m'])
    for i,f in enumerate(frames):
        s=f['state'];fluid=s['fluid'];hyd=fluid['hydraulic_points'];sampled_mechanical=fluid['mechanical_points']
        if not hyd:hyd,sampled_mechanical=reconstructed_samples(mesh,run,s['hydraulic'])
        hids=groups(hyd,len(cells));hpoints=np.asarray([q['position_m'] for q in hyd]);mpoints=np.asarray([q['position_m'] for q in mechanical]);u=np.asarray(s['mechanical']['displacement_m'])
        fields=('pressure_head_m','physical_pore_pressure_kpa','bishop_chi','effective_pore_stress_kpa','effective_saturation','water_saturation','water_content','water_equivalent_mass_per_volume')
        raw={k:np.asarray([q[k] for q in hyd]) for k in fields};raw['integration_weight_m3']=np.asarray([q['weight_m3'] for q in hyd]);raw['original_element_id']=np.asarray([q['element_index'] for q in hyd])
        constitutive={k:np.asarray([q[k] for q in sampled_mechanical]) for k in fields};constitutive.update(state_fields(s['mechanical']['material_states']))
        constitutive['integration_weight_m3']=np.asarray([q['weight_m3'] for q in mechanical]);constitutive['original_element_id']=np.asarray([q['element_index'] for q in mechanical])
        total=constitutive['effective_stress_compression_positive_kpa'].copy();total[:,:3]+=constitutive['effective_pore_stress_kpa'][:,None];constitutive['total_stress_compression_positive_kpa']=total
        cell={'original_element_id':[np.arange(len(cells))],'material_id':[np.asarray(mesh['material_ids'])],
            'cell_total_head_m':[run['case']['flow']['head_datum_m']+np.asarray(s['hydraulic']['relative_cell_head_m'])]}
        for k in fields:cell['hydraulic_volume_mean_'+k]=[average(raw[k],hyd,hids)]
        for k in ('effective_stress_compression_positive_kpa','total_stress_compression_positive_kpa','void_ratio','plastic_strain_compression_positive'):
            cell['constitutive_volume_mean_'+k]=[average(constitutive[k],mechanical,mids)]
        name=f'soil-{i:04}.vtu';point={'original_node_id':np.arange(len(xyz)),'displacement_from_seed_m':u,'displacement_from_initial_equilibrium_m':u-reference}
        meshio.write_points_cells(output/name,xyz,[(kind,cells)],point_data=point,cell_data=cell)
        checked=meshio.read(output/name);assert np.array_equal(checked.cells_dict[kind],cells)
        for k,v in point.items():assert np.array_equal(checked.point_data[k],v),k
        for k,v in cell.items():assert np.array_equal(checked.cell_data[k][0],v[0]),k
        time=str(s['hydraulic']['time_s']);ET.SubElement(collections['soil'],'DataSet',timestep=time,group='',part='0',file=name)
        for label,p,data in [('hydraulic-points',hpoints,raw),('constitutive-points',mpoints,constitutive)]:
            name=f'{label}-{i:04}.vtu';meshio.write_points_cells(output/name,p,[('vertex',np.arange(len(p))[:,None])],point_data=data);checked=meshio.read(output/name)
            assert np.array_equal(checked.points,p)
            for k,v in data.items():assert np.array_equal(checked.point_data[k],v,equal_nan=True),k
            ET.SubElement(collections[label],'DataSet',timestep=time,group='',part='0',file=name)
        balance=f['hydraulic_balance'] or {};rows.append(dict(frame=i,name=f['name'],time_s=s['hydraulic']['time_s'],water_equivalent_mass_volume_m3=fluid['water_equivalent_mass_volume_m3'],relative_force_residual=f['mechanical_audit']['relative_effective_rhs_residual'],**{k:balance.get(k,0.) for k in ('precipitation_m3_s','surface_routed_water_m3_s','storage_rate_m3_s','boundary_outflow_m3_s')}))
    for k,r in roots.items():ET.ElementTree(r).write(output/f'{k}.pvd',encoding='utf-8',xml_declaration=True)
    for k in roots:
        data=ET.parse(output/f'{k}.pvd').find('Collection');assert len(data)==len(frames)
        for row in data:assert (output/row.attrib['file']).is_file()
    with (output/'physical-time.csv').open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    receipt=dict(schema='soil-fem-one-way-export/1',passed=True,physical_frames=len(frames),raw_layouts_distinct=True,transfer_point_source='stored' if run['case']['store_transfer_integration_points'] else 'reconstructed from primary heads on original quadrature',mesh_file_sha256=sha(mesh_path),result_file_sha256=sha(result_path),solver_sha256=sha(solver),all_vtu_and_pvd_roundtrip='passed',scope='Original hydraulic and constitutive quadrature layouts remain separate. Point source is explicitly recorded. Cell means carry an explicit prefix. Soil displacement references seed and initial equilibrium separately. No nodal pressure smoothing; initial homotopy is omitted from the physical-time series.')
    with (output/'export-verification.json').open('x',encoding='utf-8') as f:json.dump(receipt,f,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('result',type=Path);p.add_argument('output',type=Path);p.add_argument('--solver',type=Path,required=True);a=p.parse_args();export(a.mesh.resolve(),a.result.resolve(),a.solver.resolve(),a.output.resolve())
