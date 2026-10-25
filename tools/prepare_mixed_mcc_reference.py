"""Common physical column for mixed MCC and the independent rate-form FEM.

Both programs read identical geometry, K0/OCR/NCL parameters and force stages.
No solved Rust state or pressure DOF initializes the independent calculation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def write(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2)


def column(n):
    xyz = [[float(x), float(y), z/n] for z in range(n+1) for y in range(2) for x in range(2)]
    tets, material = [], []
    for z in range(n):
        c = [4*z+i for i in (0, 1, 3, 2, 4, 5, 7, 6)]
        for t in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)):
            tets.append([c[i] for i in t])
            material.append(1 if 2*z < n else 2)
    edges, cells = {}, []
    for tet in tets:
        cell = tet.copy()
        for a, b in ref.EDGES:
            key = tuple(sorted((tet[a], tet[b])))
            if key not in edges:
                edges[key] = len(xyz)
                xyz.append(((np.asarray(xyz[key[0]])+xyz[key[1]])/2).tolist())
            cell.append(edges[key])
        cells.append(cell)
    xyz = np.asarray(xyz)
    regions = [dict(id=i, material_id=i, body=name, bounds_m=[[0.,0.,lo],[1.,1.,hi]], volume_m3=.5,
                    source_sha256=hashlib.sha256(b'analytic unit-footprint two-layer column').hexdigest())
               for i, name, lo, hi in ((1,'LowerClay',0.,.5),(2,'UpperClay',.5,1.))]
    return dict(schema='soil-fem-mesh/3',units='m-kN-kPa',axis='z_up',crs='LOCAL',origin_m=[0.,0.,0.],
        nodes=xyz.tolist(),tetrahedra=tets,tetrahedra10=cells,quadrature_order=3,
        material_ids=material,region_ids=material,regions=regions,face_sets={},
        node_sets=dict(all=list(range(len(xyz))),top=np.flatnonzero(xyz[:,2]==1.).tolist(),
                       bottom=np.flatnonzero(xyz[:,2]==0.).tolist(),floor=np.flatnonzero(xyz[:,2]==.5).tolist()),
        excavation_sets={})


def prepare(output):
    root = Path(__file__).resolve().parents[1]
    output.mkdir(exist_ok=False)
    ground = json.loads((root/'examples/ground-mixed-mcc.json').read_text(encoding='utf-8'))
    for material in ground['materials']:
        material['unit_weight_kn_m3'] = 20.
        material['mcc'] = {'critical_state_m':1.2,'lambda':.2,'kappa':.04,'poisson':.3,
                           'reference_pressure_kpa':1.,'specific_volume_ncl_at_reference':3.5}
    for layer in ground['layers']:
        layer.update(k0_x=1.,k0_y=1.,consolidation=dict(definition='mcc_pc_over_mean_effective',ratio=1.))
    write(output/'ground.json', ground)
    case = json.loads((root/'examples/stages-mixed-mcc.json').read_text(encoding='utf-8'))
    case['settings'].update(initial_increment=1/16,minimum_increment=1/4096,material_maximum_substep=1e-4)
    case['boundary']['constraints'] = [dict(node_set=name,direction=np.eye(3)[axis].tolist(),value_m=0.)
        for name,axis in (('all',0),('all',1),('bottom',2))]
    case['monitor_sets']={'top':'top','floor':'floor'}
    case['stages']=[dict(kind='pressure',name='load',target_pressure_kpa=130.),
        dict(kind='pressure',name='unload',target_pressure_kpa=100.),
        dict(kind='excavate',name='cut',bodies=['UpperClay'],remaining_surface_pressure_kpa=100.)]
    for n in (2,4):
        directory=output/f'n{n}';directory.mkdir()
        write(directory/'mesh.json',column(n))
        write(directory/'stages.json',case)
    write(output/'preparation.json',dict(schema='soil-fem-mixed-mcc-reference-preparation/1',
        physics='Dry isotropic K0=1, NC pc/p=1, P=100->130->100 kPa; remove upper 0.5 m, retain 100 kPa on new surface',
        reference_origin='Both start with identical input K0/NCL fields and zero displacement; no solved Rust state is used',
        limitations='Development comparison only; both spatial and integration refinements remain required',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sha256={p.relative_to(output).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in output.rglob('*.json')}))
    print('Prepared independent-input 2/4-layer mixed MCC pressure/unload/excavation columns.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path)
    prepare(parser.parse_args().output.resolve())
