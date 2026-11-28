"""New inputs using existing geometry/fields; no production feature changes."""
import copy
import json
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
import fem_reference as ref
from refine_tet10 import refine


def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def write(p, x):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(x, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def inputs(level=0, steps=4, offset=.05, beam=True):
    m, g, c = [read(ROOT / f'examples/plate6-offset-{n}.json') for n in ('mesh', 'ground', 'stages')]
    audits = []
    for _ in range(level):
        m, a = refine(m)
        audits.append(a)
    xyz = np.asarray(m['nodes'])
    m['node_sets']['bottom'] = np.flatnonzero(xyz[:, 2] == 0.).tolist()
    m['node_sets']['top'] = np.flatnonzero(xyz[:, 2] == 1.).tolist()
    cfg = c['beam_system']
    template = copy.deepcopy(cfg['plates'][0])
    topology = ref.topology(m)
    plates = []
    for i, face in enumerate(m['face_sets']['top']):
        ids = ref.face_nodes(m, topology[tuple(sorted(face))][0], face).tolist()
        if np.cross(xyz[ids[1]]-xyz[ids[0]], xyz[ids[2]]-xyz[ids[0]])[2] < 0:
            ids = [ids[j] for j in (0, 2, 1, 5, 4, 3)]
        p = copy.deepcopy(template)
        p.update(id=f'p{i}', nodes=[f'n{n}' for n in ids], attachment_offset_m=offset)
        plates.append(p)
    used = sorted({int(n[1:]) for p in plates for n in p['nodes']})
    cfg['nodes'] = [dict(id=f'n{n}', position_m=m['nodes'][n], soil_node_set=f'b02:{n}') for n in used]
    for n in used:
        m['node_sets'][f'b02:{n}'] = [n]
    cfg['plates'] = plates
    if beam:
        cfg['nodes'].append(dict(id='tip', position_m=[2., 0., 1.], match_soil_node=False))
        cfg['members'][0]['nodes'] = ['n9', 'tip']
    else:
        cfg['members'] = []
        cfg['nodal_loads'] = []
    cfg['maximum_condensed_dofs'] = 1000
    cfg['storage_budget'] = dict(maximum_nodes=500, maximum_compilation_mib=256, maximum_cached_mib=64)
    names = [p['id'] for p in plates] + [b['id'] for b in cfg['members']]
    c['stages'][0]['members'] = names
    c['stages'][-1]['members'] = names
    c['settings']['initial_increment'] = 1 / steps
    c['verification_reject_once_stages'] = []
    g['initialization']['geostatic']['verification_reject_once'] = False
    # Tighten neither physical models nor acceptance gates. Existing requested
    # solver tolerances are inherited unchanged from the old accepted example.
    return m, g, c, audits
