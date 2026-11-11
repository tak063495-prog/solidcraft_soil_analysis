"""Existing-input small A-01 specimens. No solver/default/material changes."""
import copy
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import fem_reference as fem


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, allow_nan=False, separators=(',', ':'))
        stream.write('\n')


def specimen(policy='stress_free', steps=16, quadratic=False):
    mesh = read(ROOT/'examples'/('joint-column-tet10-mesh.json' if quadratic else 'joint-column-mesh.json'))
    ground = read(ROOT/'examples/ground-joint-mc.json')
    case = read(ROOT/'examples/stages-joint-column.json')
    case['joint_system']['joints'][0]['birth_policy'] = policy
    case['settings']['initial_increment'] = 1/steps
    case['settings']['minimum_increment'] = 1/steps
    # Fixed input grid, including geostatic release. No accepted native path is
    # supplied to the independent solver.
    for stage in case['stages']:
        stage['initial_increment'] = 1/steps
    case['stages'].insert(0, dict(kind='joint_birth', name='birth', joints=['Seam'], initial_increment=1/steps))
    case['stages'] += [dict(kind='joint_death', name='death', joints=['Seam'], initial_increment=1/steps),
                       dict(kind='joint_birth', name='rebirth', joints=['Seam'], initial_increment=1/steps),
                       dict(kind='boundary', name='post_rebirth', target_factor=-.04, initial_increment=1/steps),
                       dict(kind='boundary', name='final_return', target_factor=0., initial_increment=1/steps)]
    return mesh, ground, case


def separated(policy='stress_free', steps=16, quadratic=False):
    mesh, ground, case = specimen(policy, steps, quadratic)
    n, ne = len(mesh['nodes']), len(mesh['tetrahedra'])
    base = copy.deepcopy(mesh)
    mesh['nodes'] += (np.asarray(base['nodes'])+[2., 0., 0.]).tolist()
    for key in ('tetrahedra', 'tetrahedra10'):
        if key in mesh:
            mesh[key] += (np.asarray(base[key])+n).tolist()
    for key in ('material_ids', 'region_ids'):
        mesh[key] += base[key]
    for name, nodes in base['node_sets'].items():
        if not name.startswith('Seam:'):
            mesh['node_sets'][name] += [x+n for x in nodes]
        else:
            mesh['node_sets'][name.replace('Seam:', 'Second:')] = [x+n for x in nodes]
    for name, faces in base['face_sets'].items():
        if not name.startswith('Seam:'):
            mesh['face_sets'][name] += (np.asarray(faces)+n).tolist()
        else:
            mesh['face_sets'][name.replace('Seam:', 'Second:')] = (np.asarray(faces)+n).tolist()
    for region in mesh['regions']:
        region['bounds_m'][1][0] += 2.
        region['volume_m3'] *= 2
    second = copy.deepcopy(case['joint_system']['joints'][0])
    second.update(id='Second', side_a='Second:a', side_b='Second:b')
    case['joint_system']['joints'].append(second)
    for stage in case['stages']:
        if 'joints' in stage:
            stage['joints'].append('Second')
    assert len(mesh['tetrahedra']) == 2*ne
    return mesh, ground, case


def shared(policy='stress_free', steps=16, quadratic=False):
    mesh, ground, case = specimen(policy, steps, quadratic)
    original = copy.deepcopy(case['joint_system']['joints'][0])
    case['joint_system']['joints'] = []
    for i, name in enumerate(('First', 'Second')):
        definition = copy.deepcopy(original)
        definition.update(id=name, side_a=name+':a', side_b=name+':b')
        case['joint_system']['joints'].append(definition)
        for side in ('a', 'b'):
            face = mesh['face_sets']['Seam:'+side][i]
            mesh['face_sets'][name+':'+side] = [face]
            owner, = fem.topology(mesh)[tuple(sorted(face))]
            mesh['node_sets'][name+':'+side] = fem.face_nodes(mesh, owner, face).tolist()
    for stage in case['stages']:
        if 'joints' in stage:
            stage['joints'] = ['First', 'Second']
    return mesh, ground, case


def stacked(policy='stress_free', steps=16, quadratic=False):
    """Two disjoint parallel surfaces coupled through the middle soil block."""
    from joint_tip_fixture import conformal
    from split_joint_mesh import split
    source=conformal(x=(0.,1.),y=(0.,1.),z=(0.,1/3,2/3,1.),same_body=True)
    if not quadratic:
        source.pop('tetrahedra10');source['nodes']=source['nodes'][:16]
        source['schema']='soil-fem-mesh/2';source.pop('quadrature_order')
        xyz=np.asarray(source['nodes'])
        source['node_sets']={name:np.flatnonzero(np.abs(xyz[:,2]-z)<1e-12).tolist() for name,z in (('bottom',0.),('top',1.))}
        source['node_sets']['all']=list(range(16))
    xyz=np.asarray(source['nodes']);interfaces=[]
    for name,z in (('LowerSeam',1/3),('UpperSeam',2/3)):
        faces=[]
        for face,owners in fem.topology(source).items():
            if len(owners)==2 and np.all(np.abs(xyz[list(face),2]-z)<1e-12):
                a,b=sorted(owners,key=lambda e:xyz[source['tetrahedra'][e],2].mean())
                faces.append(dict(nodes=list(face),owner_a=a,owner_b=b))
        interfaces.append(dict(id=name,faces=faces))
    mesh,_=split(source,dict(schema='soil-fem-joint-split/1',interfaces=interfaces))
    _,ground,case=specimen(policy,steps,quadratic)
    ground['materials']=ground['materials'][:1];ground['bodies']=ground['bodies'][:1]
    del ground['initialization']['geostatic']['seed']['states']['UpperClay']
    definition=case['joint_system']['joints'][0]
    case['joint_system']['joints']=[]
    for name in ('LowerSeam','UpperSeam'):
        d=copy.deepcopy(definition);d.update(id=name,side_a=name+':a',side_b=name+':b');case['joint_system']['joints'].append(d)
    for stage in case['stages']:
        if 'joints' in stage:stage['joints']=['LowerSeam','UpperSeam']
    return mesh,ground,case


def structure(kind, policy='stress_free', steps=16, quadratic=False):
    mesh, ground, case = specimen(policy, steps, quadratic)
    if kind == 'anchor':
        case['anchor_system'] = dict(schema='soil-fem-anchors/1', embedding_tolerance_m=1e-8, anchors=[dict(
            id='Tendon', young_kpa=2e6, area_m2=.001, initially_active=True,
            head=dict(position_m=[.23,.19,.8], body='UpperClay'),
            bonded_path=[dict(position_m=[.23,.19,.2], body='LowerClay'), dict(position_m=[.23,.19,.05], body='LowerClay')])])
        return mesh, ground, case
    def bound(name, position):
        selected = np.flatnonzero(np.linalg.norm(np.asarray(mesh['nodes'])-position, axis=1)<1e-12)
        # Interface coincidence is resolved by explicit node set, never by a
        # coordinate match with two possible physical sides.
        if len(selected)>1:
            selected = [n for n in selected if n in mesh['node_sets']['Seam:b']]
        assert len(selected)==1
        mesh['node_sets'][name] = [int(selected[0])]
        return dict(id=name, position_m=position, soil_node_set=name)
    cfg = dict(schema='soil-fem-structures/1', nodes=[], sections=[], members=[],
               rotation_length_scale_m=1., maximum_condensed_dofs=100,
               rotation_constraints=[], translation_constraints=[])
    if kind == 'beam':
        cfg['nodes'] = [bound('attached', [0.,0.,.5]), dict(id='fixed', position_m=[2.,0.,.5])]
        cfg['sections'] = [dict(id='S', young_kpa=2e6, poisson=.3, area_m2=.001, iy_m4=1e-5,
                                iz_m4=2e-5, torsion_j_m4=1e-5, unit_weight_kn_m3=78.5,
                                formulation=dict(kind='euler_bernoulli'))]
        cfg['members'] = [dict(id='Beam', nodes=['attached','fixed'], section_id='S', local_y=[0.,1.,0.],
                               initially_active=True, constant_line_load_kn_m=[0.,0.,-.2],
                               boundary_line_load_kn_m=[.1,.05,-.1])]
        for key in ('translation_constraints','rotation_constraints'):
            cfg[key] = [dict(node_id='fixed', direction=d.tolist(), value=0.) for d in np.eye(3)]
    elif kind == 'plate':
        cfg['nodes'] = [bound('root', [0.,0.,0.]), bound('p', [1.,0.,.5]), bound('q', [0.,1.,.5])]
        cfg['plate_sections'] = [dict(id='S', young_kpa=3e4, poisson=.2, thickness_m=.05,
                                      unit_weight_kn_m3=24., drilling_factor=1e-4)]
        cfg['plates'] = [dict(id='Plate', nodes=['root','p','q'], section_id='S', initially_active=True,
                              constant_traction_kpa=[0.,0.,-.1], boundary_traction_kpa=[.1,.05,-.2])]
        cfg['rotation_constraints'] = [dict(node_id='root', direction=d.tolist(), value=0.) for d in np.eye(3)]
    else:
        raise ValueError(kind)
    case['beam_system'] = cfg
    return mesh, ground, case


def cases(steps=16, quadratic=False):
    result = {}
    for policy in ('transfer_force', 'stress_free'):
        for name, maker in (('separated', separated), ('shared', shared)):
            result[name+'-'+policy] = maker(policy, steps, quadratic)
    for kind in ('beam', 'plate', 'anchor'):
        result['joint-'+kind] = structure(kind, steps=steps, quadratic=quadratic)
    return result
