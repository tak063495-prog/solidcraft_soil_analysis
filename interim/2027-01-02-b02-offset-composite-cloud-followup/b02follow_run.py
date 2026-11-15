"""18 native histories and 54 input-only reference histories, bounded and opt-in."""
import argparse
import copy
import json
import random
import shutil
import sys
from pathlib import Path
from b02follow_common import HERE, ROOT, OLD, legacy, write, sha, command, fingerprint, source_binding

SEED = 20270102


def permutation(n, mode, rng):
    ids = list(range(n))
    if mode == 'reverse':
        ids.reverse()
    elif mode == 'shuffle':
        rng.shuffle(ids)
    elif mode != 'default':
        raise ValueError('unknown declaration order')
    return ids


def input_case(level, scale, mode):
    if type(level) is not int or level not in (0,1) or scale not in (.5,1.,2.):
        raise ValueError('diagnostic resource/input bound: only 12/96 and the preregistered thickness series')
    m, g, c, audits = legacy['cases'].inputs(level, 4, .05*scale, True)
    original_thickness = c['beam_system']['plate_sections'][0]['thickness_m']
    for section in c['beam_system']['plate_sections']:
        section['thickness_m'] *= scale
        if section['thickness_m']/2 < .05*scale:
            raise ValueError('thickness would leave accepted offset scope')
    canonical = copy.deepcopy([m, g, c])
    rng = random.Random(SEED)
    ids = permutation(len(m['tetrahedra']), mode, rng)
    for key in ('tetrahedra', 'tetrahedra10', 'material_ids', 'region_ids'):
        m[key] = [m[key][i] for i in ids]
    faces = {}
    for key in sorted(m['face_sets']):
        faces[key] = permutation(len(m['face_sets'][key]), mode, rng)
        m['face_sets'][key] = [m['face_sets'][key][i] for i in faces[key]]
    plates = permutation(len(c['beam_system']['plates']), mode, rng)
    c['beam_system']['plates'] = [c['beam_system']['plates'][i] for i in plates]
    restored = copy.deepcopy([m, g, c])
    for key in ('tetrahedra','tetrahedra10','material_ids','region_ids'):
        restored[0][key] = [restored[0][key][ids.index(i)] for i in range(len(ids))]
    for key, order in faces.items():
        restored[0]['face_sets'][key] = [restored[0]['face_sets'][key][order.index(i)] for i in range(len(order))]
    restored[2]['beam_system']['plates'] = [restored[2]['beam_system']['plates'][plates.index(i)] for i in range(len(plates))]
    if restored != canonical:
        raise ValueError('inverse map changed physical input')
    xyz = legacy['cases'].np.asarray(m['nodes'])
    thickness = original_thickness*scale
    lengths = [max(float(legacy['cases'].np.linalg.norm(xyz[int(p['nodes'][i][1:])]-xyz[int(p['nodes'][j][1:])])) for i,j in ((0,1),(1,2),(2,0))) for p in c['beam_system']['plates']]
    mapping = dict(tetrahedra_declared_to_original=ids, faces_declared_to_original=faces,
                   plates_declared_to_original=plates, inverse_exact=True, canonical_input_sha256=fingerprint(canonical),
                   thickness_m=thickness, offset_m=.05*scale, offset_over_thickness=.5, thickness_over_L_1m=thickness,
                   thickness_over_element_max_edge=[thickness/x for x in lengths], scale=scale, mode=mode, seed=SEED)
    return m,g,c,mapping


def run(output, observer):
    output.mkdir(parents=True, exist_ok=False)
    (output/'bin').mkdir()
    binary = output/'bin'/observer.name
    shutil.copy2(observer, binary)
    write(output/'source-binding-before.json', source_binding())
    write(output/'executable-binding.json', {binary.name: sha(binary)})
    plan = []
    for level in (0,1):
        for scale in (.5,1.,2.):
            for mode in ('default','reverse','shuffle'):
                name = f'tet{12*8**level}-h{scale:g}-{mode}'
                m,g,c,mapping = input_case(level,scale,mode)
                plan.append(dict(name=name, level=level, scale=scale, mode=mode, mapping=mapping))
                folder = output/name
                folder.mkdir()
                for key,value in [('mesh',m),('ground',g),('stages',c),('order-map',mapping)]:
                    write(folder/(key+'.json'),value)
    # The complete factorial inputs and mappings exist before the first solve.
    write(output/'plan.json',dict(native_histories=18,independent_histories=54,orders=[4,5,6],cases=plan,
                                 numerical_acceptance=None,passed=False,b02_complete=False))
    for rec in plan:
        folder = output/rec['name']
        base = [folder/(n+'.json') for n in ('mesh','ground','stages')]
        command([binary,'solve',*base,folder/'native'],folder/'native.log')
        for order in (4,5,6):
            command([sys.executable,HERE/'b02follow_reference.py',folder,folder/f'gl{order}.jsonl','--order',str(order)],folder/f'gl{order}.log')
        # Existing gates only, at the original independent integration.
        shutil.copy2(folder/'gl5.jsonl',folder/'independent.jsonl')
        comparison = legacy['compare'].compare(folder)
        write(folder/'existing-gate-result.json',comparison)
        print(rec['name'], comparison['existing_field_gates_passed'],flush=True)
    after = source_binding()
    write(output/'source-binding-after.json',after)
    if json.loads((output/'source-binding-before.json').read_text()) != after:
        raise ValueError('sources changed during fresh solve')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--observer',type=Path,required=True)
    a=p.parse_args();run(a.output.resolve(),a.observer.resolve())
