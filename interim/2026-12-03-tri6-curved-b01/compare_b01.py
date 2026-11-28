"""B-02 all-field comparator extended for explicitly unattached structure nodes.

All original field gates retained; attachment checks use only defined soil nodes.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from b01_cases import read, write
from plate6_offset_reference import System
import fem_reference as ref


def metric(a, b, gate=None):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    assert a.shape == b.shape, (a.shape, b.shape)
    delta = a-b
    absolute = float(np.max(np.abs(delta), initial=0.))
    relative = float(np.linalg.norm(delta.ravel())/max(np.linalg.norm(b.ravel()), 1e-12))
    passed = None
    if gate == 'soil': passed = relative < 1e-7
    elif gate == 'kinematics': passed = absolute < 1e-9
    elif gate == 'plate': passed = bool(np.allclose(a, b, atol=1e-8, rtol=1e-7))
    elif gate == 'void': passed = absolute < 1e-8
    elif gate == 'exact': passed = bool(np.array_equal(a, b))
    return dict(maximum_absolute=absolute, relative_l2=relative, components=int(a.size), gate=gate, passed=passed)


def compare(folder):
    m, c = read(folder/'mesh.json'), read(folder/'stages.json')
    model = System(m, c['beam_system'])
    bottom = np.asarray(m['node_sets']['bottom']); nd = 3*len(m['nodes'])
    fixed = (bottom[:, None]*3+np.arange(3)).ravel(); free = np.setdiff1d(np.arange(nd), fixed)
    operators = [ref.volume(m, ei) for ei in range(len(m['tetrahedra']))]
    rows = []; maximum = {}; categories = 0
    with (folder/'native/accepted.jsonl').open() as nfile, (folder/'independent.jsonl').open() as pfile:
        for nline, pline in zip(nfile, pfile, strict=True):
            n, p = json.loads(nline), json.loads(pline); s = n['state']; fields = {}
            assert (n['accepted_index'], n['stage_index'], n['progress']) == (p['accepted_index'], p['stage_index'], p['progress'])
            def check(key, actual, expected, gate=None):
                fields[key] = metric(actual, expected, gate)
            check('soil_displacement_m', s['displacement_m'], p['displacement_m'], 'soil')
            check('soil_effective_stress_kpa', [v['effective_stress_compression_positive_kpa'] for v in s['material_states']], p['stress_kpa'], 'soil')
            strain = np.concatenate([np.einsum('qai,i->qa', q['b'], np.asarray(s['displacement_m'])[q['nodes']].ravel()) for q in operators])
            check('soil_engineering_strain_from_native_displacement', strain, p['engineering_strain'])
            check('soil_active_elements', s['active_elements'], np.ones(len(m['tetrahedra'])), 'exact')
            check('soil_material_ids', s['element_material_ids'], m['material_ids'], 'exact')
            check('soil_void_ratio', [v['void_ratio'] for v in s['material_states']], p['void_ratio'], 'void')
            check('soil_plastic_strain', [v['plastic_strain_compression_positive'] for v in s['material_states']], np.zeros((len(s['material_states']), 6)), 'exact')
            for key in ('equivalent_plastic_strain', 'plastic_work_kpa'):
                check('soil_'+key, [v[key] for v in s['material_states']], np.zeros(len(s['material_states'])), 'exact')
            assert all(v['last_return'] == 'elastic' for v in s['material_states'])
            assert all(s['active_elements'])
            for key in ('node_displacement_m', 'node_rotation_rad', 'reference_global_kinematics', 'plate_reference_global_kinematics'):
                check('structure_'+key, s['beams'][key], p['beams'][key], 'kinematics')
            for key in ('active_members', 'active_plates'):
                check('structure_'+key, s['beams'][key], p['beams'][key], 'exact')
            ni = np.asarray(n['soil_internal_kn']).copy()
            ni[np.asarray(n['interface_dofs'], int)] += n['interface_internal_kn']
            applied = np.asarray(s['applied_nodal_load_kn'])
            check('soil_internal_kn', n['soil_internal_kn'], p['soil_internal_kn'])
            check('total_internal_kn', ni, p['total_internal_kn'])
            check('applied_nodal_load_kn', applied, p['applied_nodal_load_kn'])
            reaction = np.zeros(nd); reaction[fixed] = (ni-applied)[fixed]
            check('support_reaction_by_node_kn', reaction.reshape(-1, 3), p['support_reaction_by_node_kn'])
            check('support_reaction_total_kn', reaction.reshape(-1, 3).sum(axis=0), np.asarray(p['support_reaction_by_node_kn']).sum(axis=0))
            # Observational work is evaluated with the actual native state. The
            # independently solved history p is never generated from this state.
            fb = model.assemble(s['beams'], s['boundary_factor'])[2]
            ds = np.c_[s['beams']['node_displacement_m'], s['beams']['node_rotation_rad']].ravel()
            effective = model.assemble(s['beams'], s['boundary_factor'])[3]
            attached_work = sum(float(effective[i]@np.asarray(s['displacement_m'])[node]) for i, node in enumerate(model.attach) if node is not None)
            check('load_work_kn_m', applied@np.asarray(s['displacement_m']).ravel()+fb@ds-attached_work, p['load_work_kn_m'])
            check('structural_load_work_kn_m', fb@ds, p['structural_load_work_kn_m'])
            for key in ('transition_release_load_kn',):
                if n['stage_index'] in (1, 4): check(key, s[key], p[key])
            report = n['beam_results']; pr = p['beam_results']
            for key, other in [('external_force_kn','external_force'), ('external_moment_origin_kn_m','external_moment'), ('support_force_kn','support_force'), ('support_moment_origin_kn_m','support_moment')]:
                check('structure_'+key, report[key], pr[other])
            for kind in ('members', 'plates'):
                actual = {v['id']: v for v in report.get(kind, [])}; expected = pr.get(kind, {})
                assert actual.keys() == expected.keys()
                for pid, value in expected.items():
                    def walk(a, b, path):
                        if isinstance(b, dict):
                            for k, v in b.items(): walk(a[k], v, path+'.'+k)
                        elif isinstance(b, list) and b and isinstance(b[0], dict):
                            assert len(a) == len(b)
                            for i, (aa, bb) in enumerate(zip(a, b, strict=True)): walk(aa, bb, path+f'[{i}]')
                        else: check(path, a, b, 'plate' if kind == 'plates' else None)
                    walk(actual[pid], value, kind+'.'+pid)
            for i, (nodes, section, kernel) in enumerate(model.plate_data):
                if not s['beams']['active_plates'][i]: continue
                pid = model.plates[i]['id']
                actual = next(v for v in report['plates'] if v['id'] == pid)
                expected = pr['plates'][pid]
                for qi, (aa, bb) in enumerate(zip(actual['samples'], expected['samples'], strict=True)):
                    check(f'plates.{pid}.attachment_moment_kn[{qi}]',
                          np.asarray(aa['bending_moment_kn'])-kernel.e*np.asarray(aa['membrane_force_kn_m']),
                          np.asarray(bb['bending_moment_kn'])-kernel.e*np.asarray(bb['membrane_force_kn_m']), 'plate')
                d = np.c_[s['beams']['node_displacement_m'], s['beams']['node_rotation_rad']][nodes].ravel()
                trace = 0.
                for l in ([.2,.3,.5],[.1,.8,.1],[0.,.37,.63],[.23,0.,.77],[.61,.39,0.]):
                    shape = ref.basis(np.array([[l[1],l[2],0.]]),True)[0][0][[0,1,2,4,5,6]]
                    trace = max(trace, float(np.max(np.abs(kernel.attachment_translation(l,d)-shape@d.reshape(6,6)[:,:3]))))
                for v in nodes:
                    soil = model.attach[v]
                    if soil is not None:
                        trace = max(trace, float(np.max(np.abs(np.asarray(s['beams']['node_displacement_m'])[v]-np.asarray(s['displacement_m'])[soil]))))
                fields[f'plates.{pid}.attachment_trace_m'] = dict(maximum_absolute=trace, relative_l2=0., components=15, gate='existing_attachment_1e-10_m', passed=trace<1e-10)
            residual = float(np.linalg.norm((ni-applied)[free])/max(np.linalg.norm(applied), 1.))
            for key, v in fields.items():
                group = key
                if key.startswith(('plates.', 'members.')):
                    # Retain all ID-specific rows above; compact maximum by field.
                    group = key.split('.', 1)[0]+'.'+key.split('.', 2)[2]
                    import re
                    group = re.sub(r'\[\d+\]', '[all]', group)
                old = maximum.setdefault(group, dict(maximum_absolute=0., relative_l2=0., components_compared=0, gate=v['gate'], passed=v['passed']))
                old['maximum_absolute'] = max(old['maximum_absolute'], v['maximum_absolute'])
                old['relative_l2'] = max(old['relative_l2'], v['relative_l2'])
                old['components_compared'] += v['components']
                if v['passed'] is False: old['passed'] = False
            rows.append(dict(accepted_index=n['accepted_index'], stage_index=n['stage_index'], progress=n['progress'], fields=fields,
                             native_reassembled_relative_residual=residual, independent_reassembled_relative_residual=p['free_residual_relative']))
    failures = [dict(point=r['accepted_index'], field=k) for r in rows for k, v in r['fields'].items() if v['passed'] is False]
    result = dict(accepted_points=len(rows), all_point_component_comparisons=rows, maximum_by_field=maximum,
                  existing_gate_failures=failures, existing_field_gates_passed=not failures,
                  maximum_native_reassembled_residual=max(r['native_reassembled_relative_residual'] for r in rows),
                  maximum_independent_reassembled_residual=max(r['independent_reassembled_relative_residual'] for r in rows),
                  all_field_numerical_acceptance=None, b01_complete=False,
                  note='Reaction, work and general all-field gates remain undefined; no new tolerance. Existing displacement/stress/structure/plate gates only.')
    write(folder/'comparison.json', result)
    return {k:v for k,v in result.items() if k != 'all_point_component_comparisons'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('folder', type=Path); a=p.parse_args()
    print(json.dumps(compare(a.folder), ensure_ascii=False, indent=2))
