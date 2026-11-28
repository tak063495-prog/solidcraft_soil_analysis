"""Existing direct inputs: flat Tri6 facets; no curved production element."""
import copy
import math
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
B02 = HERE.parent / '2026-12-01-tri6-offset-composite-b02'
B03 = HERE.parent / '2026-12-02-tri6-offset-kink-b03'
sys.path.insert(0, str(B02))
from cases import ROOT, inputs, read, write


def make(kind='attached', n=2, steps=4):
    m, g, c, _ = inputs(0, steps, 0., False)
    cfg = c['beam_system']
    template = copy.deepcopy(cfg['plates'][0])
    ny = 1 if kind == 'attached' else n
    vertices = []
    for j in range(ny + 1):
        for i in range(n + 1):
            if kind == 'attached':
                x = i/n
                p = [x, j/ny, 1 + math.sqrt(4-(x-.5)**2)-math.sqrt(3.75)]
            elif kind == 'roof':
                angle = math.radians(40)*i/n
                p = [25*math.sin(angle), 25*j/ny, 25*math.cos(angle)]
            else:
                raise ValueError(kind)
            vertices.append(p)
    positions = list(vertices)
    edges = {}
    triangles = []
    for j in range(ny):
        for i in range(n):
            a = j*(n+1)+i
            for tri in ((a, a+1, a+n+2), (a, a+n+2, a+n+1)):
                mids = []
                for v, w in zip(tri, (tri[1], tri[2], tri[0])):
                    edge = tuple(sorted((v, w)))
                    if edge not in edges:
                        edges[edge] = len(positions)
                        positions.append(((np.asarray(vertices[v])+vertices[w])/2).tolist())
                    mids.append(edges[edge])
                triangles.append(list(tri)+mids)
    nodes = []
    for i, p in enumerate(positions):
        node = dict(id=f'v{i}', position_m=p, match_soil_node=False)
        if kind == 'attached' and p[0] in (0., 1.):
            hits = np.flatnonzero(np.linalg.norm(np.asarray(m['nodes'])-p, axis=1) < 1e-12)
            assert len(hits) == 1
            selector = f'b01:{i}'
            m['node_sets'][selector] = [int(hits[0])]
            node['soil_node_set'] = selector
        nodes.append(node)
    cfg['nodes'] = nodes
    cfg['plates'] = [dict(copy.deepcopy(template), id=f'p{i}', nodes=[f'v{j}' for j in tri])
                     for i, tri in enumerate(triangles)]
    cfg['rotation_constraints'] = []
    cfg['translation_constraints'] = []
    cfg['maximum_condensed_dofs'] = 3000
    if kind == 'roof':
        section = cfg['plate_sections'][0]
        section.update(young_kpa=4.32e8, poisson=0., thickness_m=.25, unit_weight_kn_m3=0.)
        for plate in cfg['plates']:
            plate['boundary_traction_kpa'] = [0., 0., -90.]
        for node in nodes:
            x, y, z = node['position_m']
            tc, rc = [], []
            if x == 0.:
                tc.append([1., 0., 0.])
                rc.extend(([0., 1., 0.], [0., 0., 1.]))
            if y == 0.:
                tc.append([0., 1., 0.])
                rc.extend(([1., 0., 0.], [0., 0., 1.]))
            if y == 25.:
                tc.extend(([1., 0., 0.], [0., 0., 1.]))
                rc.append([0., 1., 0.])
            for key, directions in [('translation_constraints', tc), ('rotation_constraints', rc)]:
                for d in sorted(set(map(tuple, directions))):
                    cfg[key].append(dict(node_id=node['id'], direction=list(d), value=0.))
    names = [p['id'] for p in cfg['plates']]
    c['stages'][0]['members'] = names
    c['stages'][-1]['members'] = names
    return m, g, c, dict(kind=kind, circumferential_divisions=n, axial_divisions=ny,
                         element_geometry='straight planar facets', e_m=0.,
                         attached_nodes=sum('soil_node_set' in v for v in nodes),
                         roof_drilling_symmetry='x=0: theta_z=0; y=0: theta_x=theta_z=0',
                         production_features_added=False)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('output', type=Path)
    p.add_argument('--kind', choices=['attached', 'roof'], default='attached')
    p.add_argument('--n', type=int, default=2)
    p.add_argument('--steps', type=int, default=4)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    for key, value in zip(('mesh', 'ground', 'stages', 'geometry'), make(a.kind, a.n, a.steps)):
        write(a.output/(key+'.json'), value)
