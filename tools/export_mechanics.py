"""Export general mechanics, distinguishing raw quadrature data from averages."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import meshio
import numpy as np
from high_order_io import groups,average


def export(mesh_path,result_path,output):
    mesh=json.loads(mesh_path.read_text());result=json.loads(result_path.read_text())
    if result['schema'] not in ('soil-fem-mechanics-result/1','soil-fem-mixed-mechanics-result/1'):raise ValueError('mechanics result required')
    mixed=result['schema']=='soil-fem-mixed-mechanics-result/1'
    layout=result['integration_points'];ids=groups(layout,len(mesh['tetrahedra']))
    raw={name:np.asarray(result[name]) for name in ('stress_tension_positive_kpa','strain_engineering')}
    pressure_keys=None
    if mixed:
        pressure_keys=[(p['node_index'],p['material_id']) for p in result['pressure_dofs']]
        expected={(node,mat) for tet,mat in zip(mesh['tetrahedra'],mesh['material_ids']) for node in tet}
        if len(set(pressure_keys))!=len(pressure_keys) or set(pressure_keys)!=expected:
            raise ValueError('mixed pressure DOF inventory differs from mesh corner/material pairs')
        primary=np.array([p['compression_positive_kpa'] for p in result['pressure_dofs']]);lookup=dict(zip(pressure_keys,primary))
        gp=np.asarray(result['pressure_compression_positive_kpa'])
        if gp.shape!=(len(layout),) or not np.isfinite(primary).all() or not np.isfinite(gp).all():
            raise ValueError('invalid mixed pressure values')
        interpolated=[]
        for point in layout:
            element=point['element_index'];r=np.asarray(point['natural_coordinates']);l=np.r_[1-r.sum(),r]
            interpolated.append(l@[lookup[(node,mesh['material_ids'][element])] for node in mesh['tetrahedra'][element]])
        if not np.allclose(gp,interpolated,rtol=1e-10,atol=1e-10):raise ValueError('mixed integration pressure differs from primary interpolation')
        if not np.allclose(raw['stress_tension_positive_kpa'][:,:3].mean(axis=1),-gp,rtol=1e-10,atol=1e-10):
            raise ValueError('mixed mean stress/pressure mismatch')
        raw['pressure_compression_positive_kpa']=gp
    output.mkdir(parents=True,exist_ok=False)
    cell={'material_id':[np.asarray(mesh['material_ids'])],
          'integration_point_count':[np.asarray([len(g) for g in ids])]}
    cell.update({'volume_average_'+k:[average(v,layout,ids)] for k,v in raw.items()})
    quadratic=bool(mesh.get('tetrahedra10'));kind='tetra10' if quadratic else 'tetra'
    cells=mesh['tetrahedra10'] if quadratic else mesh['tetrahedra']
    meshio.write_points_cells(output/'mechanics.vtu',mesh['nodes'],[(kind,cells)],
        point_data={k:np.asarray(result[k]) for k in ('displacement_m','reaction_kn','applied_nodal_load_kn')},cell_data=cell)
    meshio.write_points_cells(output/'integration-points.vtu',[p['position_m'] for p in layout],
        [('vertex',np.arange(len(layout))[:,None])],point_data={**raw,
        'original_element_id':np.asarray([p['element_index'] for p in layout]),
        'integration_weight_m3':np.asarray([p['weight_m3'] for p in layout])})
    checked=meshio.read(output/'mechanics.vtu');points=meshio.read(output/'integration-points.vtu')
    assert np.array_equal(checked.cells_dict[kind],cells)
    assert np.array_equal(checked.point_data['displacement_m'],result['displacement_m'])
    for k,v in raw.items():assert np.array_equal(points.point_data[k],v)
    if mixed:
        # A shared physical corner may have two different material pressures.
        # Vertex cells retain both primary DOFs; no averaging or nodal smoothing.
        node_ids=np.array([key[0] for key in pressure_keys]);material_ids=np.array([key[1] for key in pressure_keys])
        coordinates=np.asarray(mesh['nodes'])[node_ids]
        meshio.write_points_cells(output/'pressure-dofs.vtu',coordinates,
            [('vertex',np.arange(len(primary))[:,None])],point_data=dict(original_node_id=node_ids,
                material_id=material_ids,pressure_compression_positive_kpa=primary))
        checked_pressure=meshio.read(output/'pressure-dofs.vtu')
        assert np.array_equal(checked_pressure.points,coordinates)
        for name,values in dict(original_node_id=node_ids,material_id=material_ids,pressure_compression_positive_kpa=primary).items():
            assert np.array_equal(checked_pressure.point_data[name],values)
        with (output/'pressure-dofs.csv').open('x',newline='',encoding='utf-8-sig') as f:
            writer=csv.writer(f);writer.writerow(['original_node_id','material_id','x_m','y_m','z_m','pressure_compression_positive_kpa'])
            writer.writerows([[int(node_ids[i]),int(material_ids[i]),*coordinates[i],primary[i]] for i in range(len(primary))])
        with (output/'pressure-dofs.csv').open(encoding='utf-8-sig',newline='') as f:
            rows=list(csv.DictReader(f))
        assert [(int(r['original_node_id']),int(r['material_id'])) for r in rows]==pressure_keys
        assert np.array_equal([float(r['pressure_compression_positive_kpa']) for r in rows],primary)
    record={'schema':'soil-fem-mechanics-export/1','mesh_sha256':hashlib.sha256(mesh_path.read_bytes()).hexdigest(),
        'result_sha256':hashlib.sha256(result_path.read_bytes()).hexdigest(),'cell_field_meaning':'volume_average',
        'raw_integration_point_roundtrip':'passed','displacement_roundtrip':'passed'}
    if mixed:
        record.update(formulation='Tet10/material-split P1 finite-bulk mean-pressure elasticity',
            pressure_dofs=len(pressure_keys),pressure_primary_roundtrip='passed',
            pressure_field_meaning='compression-positive mean stress; corner/material primary values and their GP interpolation, not pore pressure',
            shared_corner_material_pressures_averaged=False,
            sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file()},
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (output/'export-verification.json').write_text(json.dumps(record,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('result',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();export(a.mesh,a.result,a.output)
