"""Evidence coverage and actual facet geometry; no additional physics tolerance."""
import argparse
import json
from pathlib import Path
import numpy as np
from b01_cases import read, write
from study import exact


def run(study, probes, mutations, output):
    cases = []
    total_states = total_soil_gp = total_plate_samples = 0
    for record in read(study/'measurements.json')['cases']:
        folder = study/record['case']
        m, c = read(folder/'mesh.json'), read(folder/'stages.json')
        cfg = c['beam_system']
        positions = {v['id']: np.asarray(v['position_m']) for v in cfg['nodes']}
        normals, mid_errors = [], []
        for p in cfg['plates']:
            x = np.asarray([positions[v] for v in p['nodes']])
            n = np.cross(x[1]-x[0], x[2]-x[0])
            normals.append(n/np.linalg.norm(n))
            mid_errors.extend(np.linalg.norm(x[3:]-np.asarray([(x[i]+x[j])/2 for i, j in ((0,1),(1,2),(2,0))]), axis=1).tolist())
        geometry = dict(maximum_midpoint_distance_m=max(mid_errors),
                        minimum_normal_pair_dot=float(np.min(np.asarray(normals)@np.asarray(normals).T)),
                        all_e_zero=all(p.get('attachment_offset_m', 0.) == 0. for p in cfg['plates']),
                        structure_nodes=len(cfg['nodes']), facets=len(cfg['plates']))
        assert geometry['minimum_normal_pair_dot'] < 1., 'must have multiple genuine planes'
        if 'comparison' not in record:
            cases.append(dict(case=record['case'], geometry=geometry, completed=False))
            continue
        cmp = read(folder/'comparison.json')
        native = [json.loads(s) for s in (folder/'native/accepted.jsonl').read_text().splitlines()]
        independent = [json.loads(s) for s in (folder/'independent.jsonl').read_text().splitlines()]
        assert len(native) == len(independent) == cmp['accepted_points'] == 5*record['load_steps']
        assert [s['accepted_index'] for s in native] == list(range(1, len(native)+1))
        assert [r['accepted_index'] for r in cmp['all_point_component_comparisons']] == list(range(1, len(native)+1))
        counts = []
        for n, p, comparison in zip(native, independent, cmp['all_point_component_comparisons'], strict=True):
            assert (n['accepted_index'], n['stage_index'], n['progress']) == (p['accepted_index'], p['stage_index'], p['progress'])
            state = n['state']
            gp = len(state['material_states'])
            assert gp == len(p['stress_kpa']) == 27*len(m['tetrahedra'])
            assert len(state['beams']['node_displacement_m']) == len(cfg['nodes'])
            reports = n['beam_results'].get('plates', [])
            active = {plate['id'] for plate, flag in zip(cfg['plates'], state['beams']['active_plates'], strict=True) if flag}
            assert {v['id'] for v in reports} == active == p['beam_results']['plates'].keys()
            assert all(len(v['samples']) == 6 for v in reports)
            for report in reports:
                pid = report['id']
                for qi in range(6):
                    for field in ('membrane_force_kn_m', 'bending_moment_kn', 'shear_force_kn_m'):
                        # The comparator's recursive field names retain all IDs.
                        matching = [key for key in comparison['fields'] if key.startswith('plates.'+pid+'.') and f'[{qi}]' in key and field in key]
                        assert matching, (pid, qi, field)
            total_states += 1
            total_soil_gp += gp
            total_plate_samples += 6*len(reports)
            counts.append(dict(accepted_index=n['accepted_index'], soil_gp=gp, active_plate_ids=sorted(active), plate_samples=6*len(reports)))
        assert exact(read(folder/'native/run.json'), read(folder/'plain/run.json'))
        cases.append(dict(case=record['case'], completed=True, geometry=geometry, all_accepted_coverage=counts,
                          existing_field_gates_passed=cmp['existing_field_gates_passed']))
    restarts = read(study/'restarts.json')
    assert all(r['passed'] and r['producer_pid'] != r['resume_pid'] for r in restarts)
    scope, changed = read(probes), read(mutations)
    result = dict(passed=False, b01_complete=False, all_field_acceptance=None, cases=cases,
                  states_compared=total_states, soil_gp_states_compared=total_soil_gp,
                  active_plate_sample_states_compared=total_plate_samples,
                  restart_count=len(restarts), later_states_compared=sum(r['later_accepted_states'] for r in restarts),
                  restart_contract_passed=all(r['passed'] for r in restarts),
                  input_scope_contract_passed=scope['all_expected_contracts_passed'],
                  checkpoint_refusal_contract_passed=changed['refusal_contract_passed'],
                  evidence_coverage_passed=all(c['completed'] for c in cases))
    write(output, result)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('study', type=Path)
    p.add_argument('probes', type=Path)
    p.add_argument('mutations', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    run(a.study, a.probes, a.mutations, a.output)
