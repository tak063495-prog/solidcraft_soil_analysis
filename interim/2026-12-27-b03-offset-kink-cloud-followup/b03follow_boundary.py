"""Real Model::new boundary sweep; never replaces or adjusts its guard."""
import argparse
import copy
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np
from b03follow_math import guard, write
from b03follow_rotation import query, patch, rotations, SECTION

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def template():
    # Isolated command, so legacy generic import names never enter test discovery.
    sys.path.insert(0, str(ROOT/'interim/2026-12-02-tri6-offset-kink-b03'))
    from b03_cases import detached_fold
    mesh, _, stages, _ = detached_fold(0., .02, .02)
    cfg = stages['beam_system']
    third = copy.deepcopy(cfg['plates'][0])
    third.update(id='third', nodes=['v1', 'v0', 'third-v', 'm0-1', 'third-m02', 'third-m21'])
    cfg['plates'].append(third)
    for name, point in [('third-v', [0., -2., 1.]), ('third-m02', [0., -1., 1.]), ('third-m21', [.5, -1., 1.])]:
        cfg['nodes'].append(dict(id=name, position_m=point, match_soil_node=False))
    # No stage simulation: construction of this detached model is the observation.
    return mesh, cfg


def boundary_steps(thicknesses, ulp):
    values = {0, 1, 1024}
    # Every actual comparison is g or 2g against the CURRENT face's cutoff.
    for t in thicknesses:
        for distance_multiple in (1, 2):
            edge = 64*sys.float_info.epsilon*t/(distance_multiple*ulp)
            values.update(range(max(0, math.floor(edge)-2), math.ceil(edge)+3))
    values.update((51, 102, 153, 205, 256, 307, 409, 410, 512, 614, 819, 820))
    return sorted(values)


def run(probe, output):
    output.mkdir(parents=True, exist_ok=False)
    mesh, base_cfg = template()
    ulp = math.ulp(.02)
    tau = 64*sys.float_info.epsilon*.1
    common = ('v0', 'v1', 'm0-1')
    profiles = [('same', (.1, .1, .1)), ('mixed', (.1, .2, .05)), ('binary-exact', (.125, .125, .125))]
    all_rows, summaries = [], []
    pending_requests, pending_meta = [], []
    batch = 0

    def flush():
        nonlocal batch
        if not pending_requests:
            return
        responses = query(probe, output/'batches', str(batch).zfill(3), pending_requests)
        for meta, response in zip(pending_meta, responses, strict=True):
            meta['actual'] = response['production_model_guard']
            meta['prediction_matches'] = meta['predicted']['accepted'] == meta['actual']['accepted']
            all_rows.append(meta)
        pending_requests.clear()
        pending_meta.clear()
        batch += 1

    def add(identity, cfg, meta):
        pending_requests.append(dict(id=identity, kind='model_guard', mesh=mesh, configuration=cfg))
        pending_meta.append(dict(id=identity, **meta))
        if len(pending_requests) >= 256:
            flush()

    for profile, thicknesses in profiles:
        steps = boundary_steps(thicknesses, ulp)
        for step in steps:
            offsets = [.02, .02+step*ulp, .02+2*step*ulp]
            arms = [[0., 0., e] for e in offsets]
            for face_order in itertools.permutations(range(3)):
                prediction = guard([arms[i] for i in face_order], [thicknesses[i] for i in face_order])
                for common_order in itertools.permutations(common):
                    for private_reverse in (False, True):
                        cfg = copy.deepcopy(base_cfg)
                        nodes = {n['id']: n for n in cfg['nodes']}
                        private = [n for n in cfg['nodes'] if n['id'] not in common]
                        if private_reverse:
                            private.reverse()
                        cfg['nodes'] = [nodes[n] for n in common_order]+private
                        sections = []
                        for i, (plate, e, t) in enumerate(zip(cfg['plates'], offsets, thicknesses, strict=True)):
                            section = copy.deepcopy(cfg['plate_sections'][0])
                            section.update(id='section'+str(i), thickness_m=t)
                            sections.append(section)
                            plate.update(section_id=section['id'], attachment_offset_m=e)
                        cfg['plate_sections'] = sections
                        cfg['plates'] = [cfg['plates'][i] for i in face_order]
                        identity = f'{profile}-g{step}-p'+''.join(map(str, face_order))+'-n'+''.join(str(common.index(n)) for n in common_order)+'-r'+str(int(private_reverse))
                        add(identity, cfg, dict(profile=profile, gap_ulp=step, gap_m=step*ulp,
                            gap_over_tau=step*ulp/tau, offsets_m=offsets, thicknesses_m=thicknesses,
                            face_order=face_order, common_node_order=common_order, private_reverse=private_reverse,
                            predicted=prediction))
        # An exact integer-ULP interval enumeration of all representable steps in this bounded family.
        intervals = []
        for face_order in itertools.permutations(range(3)):
            accepted = [step for step in range(1025) if guard(
                [[0., 0., .02+i*step*ulp] for i in face_order], [thicknesses[i] for i in face_order])['accepted']]
            edge = max(accepted)
            intervals.append(dict(face_order=face_order, maximum_accepted_gap_ulp=edge,
                maximum_accepted_gap_m=edge*ulp, first_rejected_gap_ulp=edge+1,
                last_accepted_in_actual_sweep=edge in steps, first_rejected_in_actual_sweep=edge+1 in steps))
        summaries.append(dict(profile=profile, thicknesses_m=thicknesses, actual_gap_steps=steps,
            exhaustive_integer_ulp_family=[0, 1024], face_order_intervals=intervals))
    # Each face's six vertex permutations (correct midpoint permutation and signed e reversal)
    # have exactly the same z arm here. Thus all 6^3 combined local orders give the same guard.
    edge_mid = {frozenset((0, 1)): 3, frozenset((1, 2)): 4, frozenset((2, 0)): 5}
    for profile, thicknesses in profiles:
        for target_face in range(3):
            for tri in itertools.permutations(range(3)):
                cfg = copy.deepcopy(base_cfg)
                step = int(min(64*sys.float_info.epsilon*t/ulp for t in thicknesses)/2)+1
                offsets = [.02+i*step*ulp for i in range(3)]
                for i, plate in enumerate(cfg['plates']):
                    s = copy.deepcopy(cfg['plate_sections'][0]); s.update(id='s'+str(i), thickness_m=thicknesses[i])
                    cfg['plate_sections'].append(s)
                    plate.update(section_id=s['id'], attachment_offset_m=offsets[i])
                perm = list(tri)+[edge_mid[frozenset((tri[i], tri[j]))] for i, j in ((0, 1), (1, 2), (2, 0))]
                parity = 1 if tri in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1
                plate = cfg['plates'][target_face]
                plate['nodes'] = [plate['nodes'][i] for i in perm]
                plate['attachment_offset_m'] *= parity
                identity = f'local-{profile}-{target_face}-'+''.join(map(str, tri))
                add(identity, cfg, dict(profile=profile, local_permutation_face=target_face, local_vertex_order=tri,
                    gap_ulp=step, predicted=guard([[0., 0., e] for e in offsets], thicknesses)))
    # Confirm all 45 rotated patch geometries remain in the production Model construction scope.
    for cells in (1, 2, 3):
        positions, elements = patch(cells)
        for e in (0., .025, -.05):
            for label, rotation in rotations():
                cfg = copy.deepcopy(base_cfg)
                points = positions@rotation.T
                cfg['nodes'] = [dict(id='r'+str(i), position_m=p.tolist(), match_soil_node=False) for i, p in enumerate(points)]
                cfg['plate_sections'] = [dict(SECTION)]
                cfg['plates'] = [dict(id='p'+str(i), section_id=SECTION['id'], nodes=['r'+str(n) for n in nodes],
                    formulation='mindlin_tri6_offset', attachment_offset_m=e, initially_active=True) for i, nodes in enumerate(elements)]
                cfg['rotation_constraints'] = []; cfg['translation_constraints'] = []; cfg['nodal_loads'] = []
                add(f'rotation-scope-{2*cells*cells}-{e}-{label}', cfg, dict(profile='permitted-rotation-scope', predicted=dict(accepted=True)))
    flush()
    write(output/'rows.json', all_rows)
    sampled = [r for r in all_rows if 'common_node_order' in r]
    by_face = {}
    for row in sampled:
        key = (row['profile'], row['gap_ulp'], tuple(row['face_order']))
        by_face.setdefault(key, set()).add(row['actual']['accepted'])
    by_gap = {}
    for row in sampled:
        key = (row['profile'], row['gap_ulp'])
        by_gap.setdefault(key, {})[tuple(row['face_order'])] = row['actual']['accepted']
    report = dict(passed=False, b03_complete=False, diagnostic_only=True, production_scope_expanded=False,
        actual_model_constructions=len(all_rows), prediction_mismatches=[r['id'] for r in all_rows if not r['prediction_matches']],
        common_edge_node_order_sensitive_groups=sum(len(v)>1 for v in by_face.values()),
        face_order_sensitive_sampled_groups=sum(len(set(v.values()))>1 for v in by_gap.values()),
        same_coordinates_and_signed_arms_under_local_vertex_orders=True,
        arbitrary_node_list_permutation_argument='src/beams.rs: node_ids maps ID to index; positions and attachment_shifts use that same index. For detached unique IDs, a permutation only relabels slots. Private order cannot change a comparison; six shared-node orders and private reversal are actual probes, not 12! executions.',
        ulp_m=ulp, tau_m=tau, profiles=summaries,
        candidate_guard_change='Specification decision pending: symmetric all-pairs criterion with explicitly chosen thickness, or canonical representative plus diameter check. No implementation change.',
        rows_file='rows.json')
    write(output/'summary.json', report)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.probe.resolve(), a.output.resolve())
