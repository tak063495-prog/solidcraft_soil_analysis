"""Adapter for the existing checkpoint observer, strictly outside the solver."""
import copy
import json
from pathlib import Path
import numpy as np
from .compare import read, matching, validate


def configuration(original, steps, maximum=5e-6):
    case = copy.deepcopy(original)
    case['settings']['initial_increment'] = 1/steps
    case['settings']['minimum_increment'] = min(case['settings']['minimum_increment'], 1/(2*steps))
    case['verification_reject_once_stages'] = []
    for stage in case['stages']:
        stage['initial_increment'] = 1/steps
    case['joint_system']['maximum_increment_m'] = maximum
    return case


def input_subset(expected, actual, location='input'):
    """Allow serializer-added defaults, require every specified physical value."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            raise ValueError(('native input type differs', location))
        for key, value in expected.items():
            if key not in actual and (value is None or value is False or value == [] or value == {}):
                continue
            if key not in actual:
                raise ValueError(('native input missing', location+'.'+key))
            input_subset(value, actual[key], location+'.'+key)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise ValueError(('native input list differs', location))
        for i, (left, right) in enumerate(zip(expected, actual)):
            input_subset(left, right, location+f'[{i}]')
    elif expected != actual:
        raise ValueError(('native input value differs', location, expected, actual))


def convert(folder, model, steps):
    """Read actual native fields; authenticate completion, load grid and geometry."""
    folder = Path(folder); check = read(folder/'verification.json'); run = read(folder/'run.json')
    if not (check['complete'] and check['uninterrupted_run_exact_except_elapsed']
            and check['accepted_states'] == 7*steps and run['complete'] and len(run['frames']) == 7):
        raise ValueError('native path incomplete or observer verification missing')
    input_subset(model.ground, run['ground'], 'ground')
    input_subset(configuration(model.case, steps, run['case']['joint_system']['maximum_increment_m']), run['case'], 'case')
    # Reject nondefault physical additions (for example extra nodal loads)
    # that a subset comparison alone would not detect.
    from .solver import Model
    actual = Model(model.mesh, run['ground'], run['case'], pc=True)
    for field in ('k', 'seed', 'base', 'fixed', 'c', 'pressure'):
        if not np.array_equal(getattr(actual, field), getattr(model, field)):
            raise ValueError(('native physical operator differs', field))
    layouts = read(folder/'layout.json')['joint_integration_points']
    _, gp = matching(model.position, [p['position_m'] for p in layouts])
    np.testing.assert_allclose(np.array([p['normal'] for p in layouts])[gp], model.normal, atol=1e-12, rtol=0)
    np.testing.assert_allclose(np.array([p['weight_m2'] for p in layouts])[gp], model.weight, atol=1e-12, rtol=1e-11)
    stages = ['geostatic']+[s['name'] for s in model.case['stages']]
    attempts = [('geostatic', run['geostatic']['attempts'])]+[(name, [a for a in run['attempts'] if a['stage'] == name]) for name in stages[1:]]
    residuals = []; schedule_audit = []
    for name, group in attempts:
        accepted = [a for a in group if a['accepted']]
        if [a['progress_to'] for a in accepted] != [i/steps for i in range(1, steps+1)]:
            raise ValueError(('native accepted grid changed', name))
        if any(not a['committed_state_unchanged_on_rejection'] for a in group if not a['accepted']):
            raise ValueError('rejected native trial mutated committed state')
        for a in accepted:
            if a['relative_force_residual'] >= model.case['settings']['relative_force_tolerance']:
                raise ValueError('native force criterion violated')
            if any(s['true_relative_residual'] > model.case['settings']['relative_linear_tolerance'] for s in a['linear_solves']):
                raise ValueError('native true linear criterion violated')
            residuals.append(a['relative_force_residual'])
        schedule_audit.append(dict(stage=name, accepted=len(accepted), rejected=len(group)-len(accepted)))
    rows = []; fields = []
    with (folder/'accepted.jsonl').open(encoding='utf-8') as f:
        for i, line in enumerate(f):
            state = json.loads(line); phase, step = divmod(i, steps)
            if phase >= 7 or state['stage_index'] != phase:
                raise ValueError('native phase mismatch')
            progress = (step+1)/steps
            origin = 0 if phase <= 1 else model.case['stages'][phase-2]['target_factor']
            end = 0 if phase == 0 else model.case['stages'][phase-1]['target_factor']
            factor = origin+(end-origin)*progress
            if state['boundary_factor'] != factor:
                raise ValueError('native boundary factor differs from inputs')
            if any(p['equivalent_plastic_strain'] != 0 or np.any(p['plastic_strain_compression_positive']) for p in state['material_states']):
                raise ValueError('native bulk plasticity outside reference scope')
            j = state['joints']['points']
            value = {'u': np.array(state['displacement_m']),
                     'water_load': np.array(state['joint_water']['nodal_load_kn']).reshape(-1, 3),
                     'physical_pressure': np.array(state['joint_water']['physical_pressure_kpa'])[gp]}
            for key, source in [('jump', 'jump_m'), ('traction', 'traction_kpa'), ('plastic_slip', 'plastic_slip_m'), ('accumulated', 'accumulated_plastic_slip_m')]:
                value[key] = np.array([p[source] for p in j])[gp]
            value['gap'] = np.einsum('qa,qa->q', value['jump'], model.normal)
            rows.append(dict(stage=stages[phase], progress=progress, factor=factor, relative_force_residual=residuals[i]))
            fields.append(value)
    result = dict(rows=rows, fields={key: np.stack([f[key] for f in fields]) for key in fields[0]},
                  metadata=dict(steps=steps, local=None, native_contact_maximum_m=run['case']['joint_system']['maximum_increment_m'],
                                force_tolerance=model.case['settings']['relative_force_tolerance'], complete=True,
                                numerical_acceptance=False, nodes=len(model.xyz), elements=len(model.mesh['tetrahedra']),
                                contact_position=model.position.tolist(), contact_weight=model.weight.tolist(),
                                normal=model.normal.tolist(), node_position=model.xyz.tolist(), node_body=model.node_body,
                                schedule_audit=schedule_audit, separate_process_restart_verified=False))
    return validate(result)


def main():
    import argparse
    from .solver import Model
    from .compare import save, write, diagnose
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('folder', 'mesh', 'ground', 'case', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--steps', type=int, required=True); p.add_argument('--user-pc', action='store_true')
    a = p.parse_args(); model = Model(read(a.mesh), read(a.ground), read(a.case), pc=a.user_pc)
    value = convert(a.folder, model, a.steps); save(a.output, value)
    write(a.output/'diagnostics.json', diagnose(value, model))


if __name__ == '__main__':
    main()
