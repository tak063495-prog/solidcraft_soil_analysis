"""Export completed steady flow with raw GP flux distinct from cell averages."""
import argparse
import csv
import json
from pathlib import Path
import meshio
import numpy as np


def export(mesh,result,output):
    if not result['complete']:raise ValueError('only a completed saturated solution can be exported')
    output.mkdir(parents=True,exist_ok=False);active=np.asarray(result['active_elements']);xyz=np.asarray(mesh['nodes'])
    quadratic=bool(mesh.get('tetrahedra10'));cells=np.asarray(mesh['tetrahedra10'] if quadratic else mesh['tetrahedra']);kind='tetra10' if quadratic else 'tetra'
    used=np.unique(cells[active]);mapping=np.full(len(xyz),-1);mapping[used]=np.arange(len(used));head=np.asarray(result['field']['head_m'])
    assert np.all(np.asarray(result['field']['defined_nodes'])[used])
    points=result['points'];fields={k:np.asarray([p[k] for p in points]) for k in ('head_m','pressure_head_m','head_gradient','darcy_flux_m_s')}
    owners=np.asarray([p['element_index'] for p in points]);weights=np.asarray([p['weight_m3'] for p in points]);live=np.flatnonzero(active)
    means={}
    for name,values in fields.items():
        means['volume_average_'+name]=[np.asarray([np.average(values[owners==e],weights=weights[owners==e],axis=0) for e in live])]
    relative=np.asarray(result['field']['relative_head_m']);datum=np.asarray(result['field']['datum_head_m'])
    nodal={'hydraulic_head_m':head[used],'relative_head_m':relative[used],'datum_head_m':datum[used],
        'pressure_head_m':(datum[used]-xyz[used,2])+relative[used],'original_node_id':used}
    meshio.write_points_cells(output/'flow.vtu',xyz[used],[(kind,mapping[cells[active]])],point_data=nodal,cell_data={'original_element_id':[live],**means})
    meshio.write_points_cells(output/'integration-points.vtu',[p['position_m'] for p in points],[('vertex',np.arange(len(points))[:,None])],point_data={**fields,'original_element_id':owners,'integration_weight_m3':weights})
    checked=meshio.read(output/'flow.vtu');raw=meshio.read(output/'integration-points.vtu')
    for name,value in nodal.items():assert np.array_equal(checked.point_data[name],value)
    for name,value in fields.items():assert np.array_equal(raw.point_data[name],value)
    with (output/'boundary-fluxes.csv').open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(result['boundary_fluxes'][0]));writer.writeheader();writer.writerows(result['boundary_fluxes'])
    with (output/'dirichlet-nodal-outflow.csv').open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['original_node_id','reaction_outflow_m3_s'])
        writer.writerows((n,q) for n,q in enumerate(result['dirichlet_nodal_outflow_m3_s']) if q is not None)
    with (output/'water-balance.json').open('x') as stream:json.dump(result['balances'],stream,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('result',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    export(json.loads(a.mesh.read_text()),json.loads(a.result.read_text()),a.output)
